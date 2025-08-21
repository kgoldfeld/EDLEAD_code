/*
 * Bayesian Hierarchical Exchangeable (HEx) Model for 
 * Multi-Factorial Cluster Randomized Trial
 * 
 * This Stan model implements a Bayesian analysis for a 2^3 factorial design
 * within a cluster randomized trial framework. The model accounts for:
 * - Three binary treatment factors (A, B, C) and their interactions
 * - Cluster-level random effects for healthcare sites/EDs
 * - Hierarchical priors with different assumptions for main effects vs interactions
 * - Binary outcomes with logistic regression
 * 
 * Model Structure:
 * logit(P(y_i = 1)) = alpha_j[i] + tau_a*A_i + tau_b*B_i + tau_c*C_i + 
 *                     tau_ab*A_i*B_i + tau_ac*A_i*C_i + tau_bc*B_i*C_i + tau_abc*A_i*B_i*C_i
 * 
 * where:
 * - alpha_j[i] is the random effect for cluster j containing individual i
 * - tau_a, tau_b, tau_c are main effects for treatments A, B, C
 * - tau_ab, tau_ac, tau_bc are two-way interaction effects
 * - tau_abc is the three-way interaction effect
 */

data {
  
  // Sample size information
  int<lower=0> N_ED;                 // Number of clusters (emergency departments)
  int<lower=0> N;                    // Total number of patients across all clusters
  
  // Design matrix and outcomes
  matrix[N, 8] x_abc;                // Design matrix: [1, A, B, C, AB, AC, BC, ABC]
  
  // Cluster and outcome assignments
  array[N] int<lower=1,upper=N_ED> ed;     // Cluster ID for each individual
  array[N] int<lower=0,upper=1> y;         // Binary outcome (0=failure, 1=success)
  
  // Prior hyperparameters (passed from R)
  array[6] real svals;              // Vector of prior scale parameters:
                                    // [1] sigma_sigma_m: scale for simga_m prior
                                    // [2] sigma_sigma_x: scale for simga_x prior  
                                    // [3] sigma_tau_m: scale for tau_m prior
                                    // [4] sigma_tau_x: scale for tau_x prior
                                    // [5] sigma_sigma_ed: scale for sigma_alpha prior
                                    // [6] sigma_3: scale for 3-way interaction (unused)
  
}

parameters {
  
  // Non-centered parameterization for computational efficiency
  vector[8] z;                      // Standard normal draws for all coefficients
  
  // Hierarchical parameters for main effects (treatments A, B, C)
  real tau_m;                       // Common mean for main effects
  real<lower=1e-6> sigma_m;         // Standard deviation of main effects around tau_m
  
  // Hierarchical parameters for interaction effects (AB, AC, BC)
  real tau_x;                       // Common mean for 2-way interactions  
  real<lower=1e-6> sigma_x;         // Standard deviation of 2-way interactions around tau_x
  
  // Cluster-level random effects
  vector[N_ED] ed_effect;           // Random intercepts for each cluster (sigma_j)
  real<lower=1e-6> sigma_ed;        // Standard deviation of cluster effects
  
}

transformed parameters {
  
  /*
   * Non-centered parameterization transformation
   * 
   * Converts standard normal draws (z) into the actual regression coefficients.
   * This approach improves MCMC sampling efficiency by reducing correlation
   * between parameters and their scale parameters.
   * 
   * Coefficient interpretation:
   * tau[1] = tau_0    : Intercept (reference: no treatments)
   * tau[2] = tau_a    : Main effect of treatment A
   * tau[3] = tau_b    : Main effect of treatment B  
   * tau[4] = tau_c    : Main effect of treatment C
   * tau[5] = tau_ab   : Interaction effect A×B
   * tau[6] = tau_ac   : Interaction effect A×C
   * tau[7] = tau_bc   : Interaction effect B×C
   * tau[8] = tau_abc  : Three-way interaction A×B×C
   */
  
  vector[8] tau;                    // Final regression coefficients
  
  // Intercept: no hierarchical structure, just standard normal
  tau[1] = z[1];
  
  // Main effects (A, B, C): hierarchical structure with common mean tau_m
  for (i in 2:4) {
    tau[i] = sigma_m * z[i] + tau_m;
  }
  
  // Two-way interactions (AB, AC, BC): hierarchical structure with common mean tau_x
  for (i in 5:7) {
    tau[i] = sigma_x * z[i] + tau_x;
  }
  
  // Three-way interaction (ABC): fixed small scale, centered at zero
  // Note: 0.3 is a fixed scale parameter, could be made data-driven
  tau[8] = 0.3 * z[8];
  
}

model {
  
  /*
   * Prior Specifications
   * 
   * The model uses a hierarchical prior structure that allows for:
   * 1. Different prior assumptions for main effects vs interactions
   * 2. Borrowing strength across similar effect types
   * 3. Robustness through Student-t distributions for scale parameters
   */
  
  // Priors for hierarchical scale parameters (Student-t for heavy tails)
  sigma_m ~ student_t(3, 0, svals[1]);   // Scale for main effect variation
  sigma_x ~ student_t(3, 0, svals[2]);   // Scale for interaction variation
  sigma_ed ~ student_t(3, 0, svals[5]);  // Scale for cluster variation
  
  // Priors for hierarchical location parameters
  tau_m ~ normal(0, svals[3]);           // Common mean for main effects
  tau_x ~ normal(0, svals[4]);           // Common mean for interactions
  
  // Non-centered parameterization: all z parameters are standard normal
  z ~ std_normal();
  
  // Cluster random effects
  ed_effect ~ normal(0, sigma_ed);
  
  /*
   * Likelihood
   * 
   * Logistic regression with cluster random effects:
   * logit(P(y_i = 1)) = sigma_j[i] + X_i * tau
   * 
   * where sigma_j[i] is the random effect for the cluster containing individual i
   * and X_i * tau represents the linear combination of treatment effects
   */
  y ~ bernoulli_logit(ed_effect[ed] + x_abc * tau);
  
}

generated quantities {
  
  /*
   * Treatment Effect Contrasts
   * 
   * Computes log-odds ratios for each treatment combination relative to control.
   * These quantities facilitate interpretation and comparison of treatment effects.
   * 
   * Each lOR represents the log-odds ratio comparing a specific treatment
   * combination to the control condition (no treatments).
   */
  
  array[7] real lOR;                // Log-odds ratios for treatment combinations
  
  // Single treatment effects (vs. control)
  lOR[1] = tau[2];                                            // A vs None
  lOR[2] = tau[3];                                            // B vs None  
  lOR[3] = tau[4];                                            // C vs None
  
  // Two-way treatment combinations (vs. control)
  lOR[4] = tau[2] + tau[3] + tau[5];                          // AB vs None
  lOR[5] = tau[2] + tau[4] + tau[6];                          // AC vs None
  lOR[6] = tau[3] + tau[4] + tau[7];                          // BC vs None
  
  // Three-way treatment combination (vs. control)
  lOR[7] = tau[2] + tau[3] + tau[4] + tau[5] + tau[6] + tau[7] + tau[8];  // ABC vs None
  
  /*
   * Additional quantities that could be computed:
   * 
   * // Convert to odds ratios
   * array[7] real OR;
   * for (i in 1:7) OR[i] = exp(lOR[i]);
   * 
   * // Posterior predictive checks
   * array[N] int y_rep;
   * for (i in 1:N) {
   *   y_rep[i] = bernoulli_logit_rng(ed_effect[ed[i]] + x_abc[i] * tau);
   * }
   * 
   * // Treatment effect comparisons (e.g., A vs B)
   * real lOR_A_vs_B = lOR[1] - lOR[2];
   */
  
}

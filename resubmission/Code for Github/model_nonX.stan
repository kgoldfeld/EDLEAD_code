/*
 * Simple (Non-exchangeable) Model for Multi-Factorial Cluster Randomized Trial
 * 
 * This Stan model provides a simpler alternative to the hierarchical exchange (HEx)
 * model for analyzing a 2^3 factorial cluster randomized trial. Unlike the HEx model,
 * this approach treats all treatment effects as exchangeable with identical priors,
 * without distinguishing between main effects and interactions.
 * 
 * Key Differences from HEx Model:
 * - No hierarchical structure for treatment effects
 * - All coefficients (main effects + interactions) have identical priors
 * - Simpler parameterization with fewer hyperparameters
 * - No non-centered parameterization needed
 * 
 * Model Structure:
 * logit(P(y_i = 1)) = α_j[i] + τ_a*A_i + τ_b*B_i + τ_c*C_i + 
 *                     τ_ab*A_i*B_i + τ_ac*A_i*C_i + τ_bc*B_i*C_i + τ_abc*A_i*B_i*C_i
 * 
 * 
 * Prior Structure:
 * τ_k ~ Normal(0, σ_τ) for k = a,...,abc, σ_τ is user-specified
 * α_j ~ Normal(0, σ_α) for j = 1,...,J , σ_α is user-specified
 * 
 * Use Cases:
 * - When you want to treat all treatment effects as equally likely a priori
 * - Computational simplicity and faster sampling
 * - Sensitivity analysis comparing hierarchical vs non-hierarchical approaches
 * - When domain knowledge doesn't suggest different priors for main effects vs interactions
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
                                    // [1-4] unused (for compatibility with HEx model)
                                    // [5] sigma_sigma_ed: scale for σ_α prior (cluster effects)
                                    // [6] sigma_tau: scale for treatment effect priors
  
}

parameters {
  
  /*
   * Direct Parameterization (vs Non-Centered in HEx Model)
   * 
   * This model uses direct parameterization which is simpler but may be
   * less computationally efficient than the non-centered approach in the HEx model.
   * However, for this relatively simple structure, the efficiency difference
   * is likely minimal.
   */
  
  // Treatment effect coefficients (all treated identically)
  vector[8] tau;                    // [τ₀, τ₁, τ₂, τ₃, τ₄, τ₅, τ₆, τ₇]
                                    // τ₀: Intercept
                                    // τ₁-τ₃: Main effects (A, B, C)
                                    // τ₄-τ₆: Two-way interactions (AB, AC, BC)
                                    // τ₇: Three-way interaction (ABC)
  
  // Cluster-level random effects (identical to HEx model)
  vector[N_ED] ed_effect;           // Random intercepts for each cluster (α_j)
  real<lower=0> sigma_ed;           // Standard deviation of cluster effects
  
}

model {
  
  /*
   * Prior Specifications
   * 
   * This model uses a much simpler prior structure compared to the HEx model:
   * - All treatment effects have identical Normal(0, σ_τ) priors
   * - No distinction between main effects and interactions
   * - No hierarchical borrowing of strength across effect types
   */
  
  // Cluster random effects (same as HEx model)
  ed_effect ~ normal(0, sigma_ed);
  sigma_ed ~ student_t(3, 0, svals[5]);  // Scale for cluster variation
  
  // Treatment effects: all coefficients treated identically
  // Note: This includes the intercept tau[1], which gets Normal(0, σ_τ) prior
  tau ~ normal(0, svals[6]);             // Common prior for all treatment effects
  
  /*
   * Likelihood (identical to HEx model)
   * 
   * Logistic regression with cluster random effects:
   * logit(P(y_i = 1)) = α_j[i] + X_i * τ
   */
  y ~ bernoulli_logit(ed_effect[ed] + x_abc * tau);
  
}

generated quantities {
  
  /*
   * Treatment Effect Contrasts (identical to HEx model)
   * 
   * Computes log-odds ratios for each treatment combination relative to control.
   * The interpretation is identical to the HEx model, but the underlying
   * coefficients come from different prior structures.
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
   * COMPARISON WITH HEx MODEL:
   * 
   * Advantages of Simple Model:
   * - Fewer hyperparameters to specify
   * - Simpler to understand and explain
   * - Faster computation
   * - Less risk of prior specification issues
   * 
   * Advantages of HEx Model:
   * - Can incorporate domain knowledge about effect hierarchies
   * - Borrows strength across similar effect types
   * - More flexible prior structure
   * - Better computational properties (non-centered parameterization)
   * 
   * PRACTICAL CONSIDERATIONS:
   * 
   * Choose Simple Model when:
   * - You have no strong prior beliefs about effect hierarchies
   * - You want computational simplicity
   * - You're doing sensitivity analysis
   * - Sample size is adequate for estimating all effects
   *

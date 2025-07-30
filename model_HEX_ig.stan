/*
 * HEx Inverse Gamma Prior Model for Multi-Factorial Cluster Randomized Trial
 * 
 * This Stan model implements an alternative hierarchical approach to the HEx model,
 * using inverse gamma priors on variance parameters instead of Student-t priors
 * on standard deviation parameters. This represents a different approach to
 * specifying prior beliefs about the scale of treatment effect heterogeneity.
 * 
 * Key Differences from Standard HEx Model:
 * - Inverse gamma priors on variance parameters (σ²) instead of Student-t on SD (σ)
 * - Non-centered parameterization with variance-based transformation
 * - Direct variance estimation with derived standard deviations
 * 
 * Model Structure (identical to HEx):
 * logit(P(y_i = 1)) = α_j[i] + τ_a*A_i + τ_b*B_i + τ_c*C_i + 
 *                     τ_ab*A_i*B_i + τ_ac*A_i*C_i + τ_bc*B_i*C_i + τ_abc*A_i*B_i*C_i
 * 
 * Prior Structure:
 * σ²_main ~ InverseGamma(α₁, β₁)     [Variance for main effects]
 * σ²_int ~ InverseGamma(α₂, β₂)      [Variance for interactions]  
 * σ²_cluster ~ InverseGamma(α₃, β₃)  [Variance for cluster effects]
 * 
 * Use Cases:
 * - When you have conjugate prior beliefs about variance parameters
 * - Sensitivity analysis with different prior families
 * - When inverse gamma priors better reflect domain knowledge
 * - Comparison with Student-t based approaches
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
  
  // Prior hyperparameters for inverse gamma distributions
  array[6] real svals;              // Vector of inverse gamma parameters:
                                    // [1] α₁: shape parameter for main effect variance
                                    // [2] β₁: scale parameter for main effect variance
                                    // [3] σ_τm: scale for main effect location prior
                                    // [4] σ_τx: scale for interaction location prior
                                    // [5] α₂: shape parameter for cluster variance
                                    // [6] β₂: scale parameter for cluster variance
                                    // Note: Uses same shape/scale for main & interaction variances
  
}

parameters {
  
  /*
   * Non-Centered Parameterization with Variance Parameters
   * 
   * This approach parameterizes the model in terms of variances (σ²) rather
   * than standard deviations (σ), then derives SDs in transformed parameters.
   * This can have different computational properties and prior implications.
   */
  
  // Standard normal draws for non-centered parameterization
  vector[8] z;                      // Standard normal draws for all coefficients
  
  // Hierarchical location parameters (identical to HEx model)
  real tau_m;                       // Common mean for main effects
  real tau_x;                       // Common mean for interaction effects
  
  // Variance parameters (key difference from HEx model)
  real<lower=0> var_m;              // Variance of main effects around tau_m
  real<lower=0> var_x;              // Variance of interactions around tau_x
  real<lower=0> var_ed;             // Variance of cluster random effects
  
  // Cluster-level random effects
  vector[N_ED] ed_effect;           // Random intercepts for each cluster (α_j)
  
}

transformed parameters {
  
  /*
   * Variance-to-Standard Deviation Transformation
   * 
   * Convert variance parameters to standard deviations for use in
   * the non-centered parameterization. This ensures computational
   * compatibility with the coefficient transformations.
   */
  
  // Derived standard deviations from variance parameters
  real<lower=0> sigma_m = sqrt(var_m);    // SD for main effects
  real<lower=0> sigma_x = sqrt(var_x);    // SD for interactions
  real<lower=0> sigma_ed = sqrt(var_ed);  // SD for cluster effects
  
  /*
   * Non-Centered Parameterization Transformation
   * 
   * Identical to HEx model: converts standard normal draws into
   * properly scaled regression coefficients using hierarchical structure.
   */
  
  vector[8] tau;                    // Final regression coefficients
  
  // Intercept: just the standard normal draw
  tau[1] = z[1];
  
  // Main effects (A, B, C): hierarchical structure with variance-derived SD
  for (i in 2:4) {
    tau[i] = sigma_m * z[i] + tau_m;
  }
  
  // Two-way interactions (AB, AC, BC): hierarchical structure
  for (i in 5:7) {
    tau[i] = sigma_x * z[i] + tau_x;
  }
  
  // Three-way interaction (ABC): fixed scale (same as HEx model)
  tau[8] = 0.3 * z[8];
  
}

model {
  
  /*
   * Inverse Gamma Prior Specifications
   * 
   * Key difference from HEx model: uses inverse gamma priors on variance
   * parameters rather than Student-t priors on standard deviations.
   * 
   * Inverse Gamma Properties:
   * - Conjugate for normal likelihood variance parameters
   * - Right-skewed distribution (mass concentrated near zero)
   * - Mean = β/(α-1) for α > 1
   * - Mode = β/(α+1)
   * - Can be more/less informative than Student-t depending on parameters
   */
  
  // Inverse gamma priors on variance parameters
  // Note: Using target += for efficiency with log probability density
  target += inv_gamma_lpdf(var_m | svals[1], svals[2]);   // Main effect variance
  target += inv_gamma_lpdf(var_x | svals[1], svals[2]);   // Interaction variance  
  target += inv_gamma_lpdf(var_ed | svals[5], svals[6]);  // Cluster variance
  
  // Normal priors on hierarchical location parameters (same as HEx)
  tau_m ~ normal(0, svals[3]);      // Common mean for main effects
  tau_x ~ normal(0, svals[4]);      // Common mean for interactions
  
  // Standard components (identical to HEx model)
  ed_effect ~ normal(0, sigma_ed);  // Cluster random effects
  z ~ std_normal();                 // Non-centered parameterization
  
  /*
   * Likelihood (identical to HEx model)
   */
  y ~ bernoulli_logit(ed_effect[ed] + x_abc * tau);
  
}

generated quantities {
  
  /*
   * Treatment Effect Contrasts (identical to other models)
   * 
   * Despite different prior specification, the interpretation of
   * treatment effects remains the same across all model variants.
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
   * COMPARISON WITH OTHER MODELS:
   * 
   * vs HEx Model (Student-t priors):
   * --------------------------------
   * Inverse Gamma Model:
   * + Conjugate priors may improve computational efficiency
   * + Natural for variance parameters
   * + Can encode strong beliefs about small variances
   * - Less robust to prior misspecification
   * - Right-skewed may not reflect uncertainty well
   * 
   * HEx Model:
   * + More robust, heavy-tailed priors
   * + Student-t more commonly used in practice
   * + Better for weakly informative priors
   * - Non-conjugate (though rarely matters with MCMC)
   * 
   * vs Simple Model:
   * ----------------
   * Inverse Gamma Model:
   * + Hierarchical structure borrows strength
   * + Separate priors for main effects vs interactions
   * + Better regularization with many parameters
   * - More complex, more hyperparameters
   * - Requires careful prior elicitation
   * 
   * PRACTICAL CONSIDERATIONS:
   * 
   * Choose Inverse Gamma Model when:
   * - You have strong prior beliefs about variance magnitudes
   * - You want conjugate prior structure
   * - Your domain suggests right-skewed variance priors
   * - You're comparing different prior families
   * 
   * Prior Elicitation for Inverse Gamma:
   * - Shape α: controls concentration around mode
   * - Scale β: controls location of distribution
   * - Higher α = more concentrated around mode
   * - Consider prior predictive checks to validate choices
   * 
   * SENSITIVITY ANALYSIS:
   * Compare results across Student-t, inverse gamma, and simple models
   * to assess robustness of conclusions to prior specification.
   */
  
}

/*
 * Prior Predictive Model for Multi-Factorial Cluster Randomized Trial
 * 
 * This Stan model generates samples from the prior predictive distribution
 * for a 2^3 factorial cluster randomized trial. It's used for prior checking
 * to ensure that prior specifications lead to reasonable predicted outcomes
 * before fitting the model to real data.
 * 
 * Purpose:
 * - Validate prior assumptions by examining implied outcome distributions
 * - Check if priors allow for clinically plausible treatment effects
 * - Identify potential issues with prior specifications
 * - Guide prior elicitation through visualization of implied predictions
 * 
 * Workflow:
 * 1. Sample all model parameters from their prior distributions
 * 2. Generate predicted outcomes (y_rep) for each observation
 * 3. Analyze distributions of predicted probabilities by treatment arm
 * 
 * This model mirrors the main analysis model but operates in "prior-only" mode.
 */

data {
  
  // Sample size information (same as main model)
  int<lower=0> N_ED;                 // Number of clusters (emergency departments)
  int<lower=0> N;                    // Total number of patients
  
  // Design structure (same as main model)
  matrix[N, 8] x_abc;                // Design matrix: [1, A, B, C, AB, AC, BC, ABC]
  array[N] int<lower=1,upper=N_ED> ed;     // Cluster assignment for each individual
  array[N] int<lower=0,upper=1> y;         // Observed outcomes (unused in prior predictive)
  
  // Prior hyperparameters (explicitly named for clarity)
  real sigma_tau_m;                  // Prior SD for main effect location parameters
  real sigma_tau_x;                  // Prior SD for interaction location parameters  
  real sigma_sigma_m;                // Prior scale for main effect variation
  real sigma_sigma_x;                // Prior scale for interaction variation
  real sigma_3;                      // Prior scale for three-way interaction
  real sigma_sigma_ed;               // Prior scale for cluster effect variation
  
}

generated quantities {
  
  /*
   * STEP 1: Sample Hierarchical Location Parameters
   * 
   * These represent the "typical" or "average" treatment effects
   * around which individual effects vary.
   */
  
  // Common mean for main effects (A, B, C)
  real tau_m = normal_rng(0, sigma_tau_m);
  
  // Common mean for two-way interactions (AB, AC, BC)  
  real tau_x = normal_rng(0, sigma_tau_x);
  
  /*
   * STEP 2: Sample Hierarchical Scale Parameters
   * 
   * These control how much individual effects vary around their common means.
   * Using while loops ensures positive values since Student-t can be negative.
   * This is equivalent to half-Student-t priors but implemented via rejection sampling.
   */
  
  // Scale parameter for main effect variation
  real sigma_m = -1;
  while (sigma_m < 0) {
    sigma_m = student_t_rng(3, 0, sigma_sigma_m);
  }
  
  // Scale parameter for interaction effect variation
  real sigma_x = -1;
  while (sigma_x < 0) {
    sigma_x = student_t_rng(3, 0, sigma_sigma_x);
  }
  
  // Scale parameter for cluster random effects
  real sigma_ed = -1;
  while (sigma_ed < 0) {
    sigma_ed = student_t_rng(3, 0, sigma_sigma_ed);
  }
  
  /*
   * STEP 3: Sample Cluster Random Effects
   * 
   * Each cluster gets its own random intercept representing
   * cluster-specific factors affecting baseline risk.
   */
  
  array[N_ED] real ed_effect;
  for (i in 1:N_ED) {
    ed_effect[i] = normal_rng(0, sigma_ed);
  }
  
  /*
   * STEP 4: Sample Standard Normal Draws for Non-Centered Parameterization
   * 
   * These will be transformed into the actual regression coefficients
   * using the hierarchical structure.
   */
  
  vector[8] z;
  for (i in 1:8) {
    z[i] = normal_rng(0, 1);
  }
  
  /*
   * STEP 5: Transform to Final Regression Coefficients
   * 
   * This mirrors the transformed parameters block in the main model,
   * converting standard normal draws into properly scaled coefficients.
   */
  
  vector[8] tau;                    // Final regression coefficients
  
  // Intercept: just the standard normal draw
  tau[1] = z[1];
  
  // Main effects (A, B, C): hierarchical structure
  for (i in 2:4) {
    tau[i] = sigma_m * z[i] + tau_m;
  }
  
  // Two-way interactions (AB, AC, BC): hierarchical structure  
  for (i in 5:7) {
    tau[i] = sigma_x * z[i] + tau_x;
  }
  
  // Three-way interaction (ABC): directly scaled (no hierarchical mean)
  tau[8] = sigma_3 * z[8];
  
  /*
   * STEP 6: Compute Treatment Effect Contrasts
   * 
   * Calculate log-odds ratios for each treatment combination
   * relative to the control condition. These are the key quantities
   * of interest for understanding treatment effects.
   */
  
  array[7] real lOR;                // Log-odds ratios vs control
  
  // Single treatments vs control
  lOR[1] = tau[2];                                            // A vs None
  lOR[2] = tau[3];                                            // B vs None
  lOR[3] = tau[4];                                            // C vs None
  
  // Two-way combinations vs control
  lOR[4] = tau[2] + tau[3] + tau[5];                          // AB vs None
  lOR[5] = tau[2] + tau[4] + tau[6];                          // AC vs None
  lOR[6] = tau[3] + tau[4] + tau[7];                          // BC vs None
  
  // Three-way combination vs control
  lOR[7] = tau[2] + tau[3] + tau[4] + tau[5] + tau[6] + tau[7] + tau[8];  // ABC vs None
  
  /*
   * STEP 7: Generate Prior Predictive Outcomes
   * 
   * For each observation, compute the linear predictor and generate
   * a binary outcome from the implied Bernoulli distribution.
   * These y_rep values show what outcomes the prior expects to see.
   */
  
  array[N] int y_rep;               // Prior predictive outcomes
  array[N] real linpred;            // Linear predictors (for debugging/analysis)
  
  for (n in 1:N) {
    // Linear predictor: cluster effect + treatment effects
    linpred[n] = ed_effect[ed[n]] + dot_product(row(x_abc, n), tau);
    
    // Generate binary outcome from logistic model
    y_rep[n] = bernoulli_logit_rng(linpred[n]);
  }
  
  /*
   * INTERPRETATION OF OUTPUTS:
   * 
   * y_rep: Shows what outcomes we expect to see under the prior
   *        - Can compute treatment arm-specific success rates
   *        - Should span reasonable clinical ranges
   *        - Extreme values may indicate problematic priors
   * 
   * lOR: Shows the range of treatment effects implied by priors
   *      - Values near 0 suggest small effects
   *      - Large absolute values suggest strong effects
   *      - Distribution should match clinical expectations
   * 
   * linpred: Linear predictors on log-odds scale
   *          - Can help diagnose numerical issues
   *          - Extreme values may cause computational problems
   * 
   * PRIOR CHECKING WORKFLOW:
   * 1. Run this model with candidate prior specifications
   * 2. Examine distributions of y_rep by treatment arm
   * 3. Check if implied success rates are clinically reasonable
   * 4. Adjust priors if predictions are unrealistic
   * 5. Iterate until satisfied with prior implications
   */
  
}

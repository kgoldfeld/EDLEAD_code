#' ---
#' title: "Bayesian Power Analysis for Sample Size Determination"
#' author: "Keith Goldfeld"
#' date: "`r Sys.Date()`"
#' description: "Monte Carlo simulation study to determine required sample size
#'              for adequate precision in Bayesian analysis of factorial cluster
#'              randomized trial based on posterior standard deviation criteria"
#' ---

# =============================================================================
# REQUIRED LIBRARIES
# =============================================================================

library(simstudy)     # Data simulation
library(data.table)   # Efficient data manipulation
library(cmdstanr)     # Bayesian modeling with Stan
library(posterior)    # Working with posterior draws
library(parallel)     # Parallel computing for Monte Carlo simulation
library(ggplot2)      # Data visualization

#' Bayesian Model Fitting and Posterior Summarization
#' 
#' Fits a Bayesian model to simulated trial data and extracts key posterior
#' statistics for power/precision analysis.
#' 
#' @param generated_data A data.table containing simulated trial data from s_generate()
#' @param mod A compiled Stan model object
#' @param argsvec Named vector containing hyperparameter values for the model
#' 
#' @return A data.table with posterior summaries for each parameter:
#'   \itemize{
#'     \item{variable: parameter name (lOR, tau, sigma, delta)}
#'     \item{p.025, p.25, p.50, p.75, p.95, p.975: posterior quantiles}
#'     \item{sd: posterior standard deviation (key metric for precision)}
#'     \item{var: posterior variance}
#'     \item{p_less_zero: probability parameter is negative}
#'     \item{p_meaningful: probability of clinically meaningful effect (< -0.223)}
#'   }
#' 
#' @details The function:
#'   1. Converts simulated data to Stan format
#'   2. Fits Bayesian model with robust MCMC settings
#'   3. Extracts posterior draws for treatment effects and model parameters
#'   4. Computes posterior summaries including precision metrics
#'   
#'   MCMC Settings optimized for complex hierarchical models:
#'   - 4 chains with 500 warmup + 2500 sampling iterations
#'   - High adapt_delta (0.98) and max_treedepth (20) for difficult geometry
#'   - Parallel chain execution for efficiency

s_bayes <- function(generated_data, mod, argsvec) {
  
  #' Convert Data.Table to Stan Data List
  #' 
  #' Internal helper function that transforms simulated trial data
  #' into the list format required by Stan models.
  
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                               # Total number of observations 
    x_abc <- model.matrix(~r1*r2*r3, data = dx) # Full factorial design matrix
    y <- dx[, y]                                # Binary outcomes
    N_ED <- dx[, length(unique(ed))]            # Number of clusters
    ed <- dx[, ed]                              # Cluster identifiers
    svals <- c(s1, s2, s3, s4, s5, s6)          # Prior hyperparameters
    
    list(N_ED = N_ED, N = N, x_abc = x_abc, ed = ed, y = y, svals = svals)
  }
  
  # Import hyperparameters into local environment
  list2env(as.list(argsvec), envir = environment())
  
  # Fit Bayesian model with robust MCMC settings
  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,                    # Suppress progress output for batch processing
    chains = 4L,                    # Multiple chains for convergence assessment
    parallel_chains = 4L,           # Parallel execution for efficiency
    iter_warmup = 500,              # Warmup iterations for adaptation
    iter_sampling = 2500,           # Sampling iterations per chain
    adapt_delta = 0.98,             # High acceptance rate for difficult geometry
    max_treedepth = 20,             # Allow deep trees for complex posteriors
    show_messages = FALSE           # Suppress Stan messages for batch processing
  )
  
  # Convert to posterior draws format
  posterior <- as_draws_array(fit$draws())
  
  # Extract and summarize key parameters
  # Focus on treatment effects (lOR), coefficients (tau), and scale parameters (sigma)
  sumbayes <- as.data.table(posterior)[
    substr(variable, 1, 3) == "lOR" |          # Log-odds ratios (treatment effects)
      substr(variable, 1, 3) == "tau" |        # Regression coefficients  
      substr(variable, 1, 5) == "sigma" |      # Scale parameters
      substr(variable, 1, 5) == "delta",       # Additional parameters (if present)
    .(
      p.025 = quantile(value, 0.025),          # 95% credible interval bounds
      p.25 = quantile(value, 0.25),            # IQR bounds
      p.50 = quantile(value, 0.50),            # Posterior median
      p.75 = quantile(value, 0.75),            # IQR bounds
      p.95 = quantile(value, 0.95),            # 90% credible interval
      p.975 = quantile(value, 0.975),          # 95% credible interval bounds
      sd = sd(value),                          # PRECISION METRIC: posterior SD
      var = var(value),                        # Posterior variance
      p_less_zero = mean(value < 0),           # Probability of negative effect
      p_meaningful = mean(value < -0.223)      # Probability of meaningful effect
      # (log(0.8) ≈ -0.223 corresponds to OR < 0.8)
    ), keyby = variable]
  
  return(sumbayes)  # Return posterior summary table
}

#' Single Replication of Power Analysis
#' 
#' Performs one complete replication of the power analysis: generates data,
#' fits model, and extracts results.
#' 
#' @param argsvec Named vector of simulation parameters including:
#'   \itemize{
#'     \item{n_ed: number of clusters}
#'     \item{icc: intracluster correlation}
#'     \item{grp1-grp4: cluster size parameters}
#'     \item{t_0, t_a, t_b, t_c: effect parameters}
#'     \item{s1-s6: prior hyperparameters}
#'   }
#' @param mod Compiled Stan model object
#' 
#' @return List containing:
#'   \itemize{
#'     \item{args: input simulation parameters}
#'     \item{model_res: posterior summary table from s_bayes()}
#'   }
#' 
#' @details This function represents one Monte Carlo iteration:
#'   1. Generate simulated dataset under null hypothesis
#'   2. Fit Bayesian model to simulated data  
#'   3. Extract posterior summaries
#'   4. Return results with parameter tracking

s_replicate <- function(argsvec, mod) {
  
  # Generate simulated dataset
  # Note: Requires s_define() and s_generate() functions from previous code
  list_of_defs <- s_define()                     # Create data definitions
  generated_data <- s_generate(list_of_defs, argsvec)  # Generate one dataset
  
  # Fit Bayesian model and extract posterior summaries
  model_bayes_1 <- s_bayes(generated_data, mod, argsvec)
  
  # Package results with input parameters for tracking
  summary_stats <- c(
    list(args = argsvec,              # Input parameters
         model_res = model_bayes_1)   # Posterior summaries
  )
  
  return(summary_stats)
}

#' Create Parameter Scenarios for Power Analysis
#' 
#' Generates all combinations of specified parameter values.
#' Identical to previous definition but included for completeness.

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

# =============================================================================
# POWER ANALYSIS PARAMETERS
# =============================================================================

# Sample size range for evaluation
# Testing 6 different cluster counts: 48, 56, 64, 72, 80, 88 clusters
n_ed <- 8 * c(6, 7, 8, 9, 10, 11)

# Fixed study design parameters
icc <- 0.015                         # Intracluster correlation coefficient
n_quarters <- 2                      # Time periods per cluster

# Base cluster sizes (before scaling)
grp1 <- 60                          # Small clusters
grp2 <- 90                          # Medium-small clusters  
grp3 <- 190                         # Medium-large clusters
grp4 <- 300                         # Large clusters

# Cluster size scaling based on expected patient volume
# Testing scenario with 40 patients/month average
mean <- 40                          # Expected patients per month per cluster
# Alternative values: 25, 30, 35, 40
ratio <- mean/40                    # Scaling factor from base scenario

# Apply scaling to all cluster size groups
grp1 <- grp1 * ratio
grp2 <- grp2 * ratio  
grp3 <- grp3 * ratio
grp4 <- grp4 * ratio

# Treatment effect parameters (NULL HYPOTHESIS scenario)
# All effects set to zero to evaluate Type I error control and precision
t_0 <- -0.40                        # Baseline log-odds (approximately 40% control rate)
t_a <- 0                           # Main effect A (null)
t_b <- 0                           # Main effect B (null)  
t_c <- 0                           # Main effect C (null)
x_ab <- 0                          # Interaction A×B (null)
x_ac <- 0                          # Interaction A×C (null)
x_bc <- 0                          # Interaction B×C (null)
x_abc <- 0                         # Three-way interaction A×B×C (null)

# Prior hyperparameter scenario
# These values determine the informativeness of priors
sds <- c(0.3, 0.3, 0.4, 0.4, 0.5, 0.4)  # Scenario 1 hyperparameters
# [s1, s2, s3, s4, s5, s6]

# Generate all parameter combinations
scenarios <- scenario_list(
  n_ed = n_ed, icc = icc, n_quarters = n_quarters,
  grp1 = grp1, grp2 = grp2, grp3 = grp3, grp4 = grp4,
  t_0 = t_0, t_a = t_a, t_b = t_b, t_c = t_c, 
  x_ab = x_ab, x_ac = x_ac, x_bc = x_bc, x_abc = x_abc,
  s1 = sds[1], s2 = sds[2], s3 = sds[3], s4 = sds[4], s5 = sds[5], s6 = sds[6]
)

# Replicate each scenario for Monte Carlo simulation
# 1000 replications per parameter combination for stable estimates
# but be forewarned that this was actually run on an HPC. If you
# want code to set that up, pleast contact me.
scenarios <- rep(scenarios, each = 1000)

# =============================================================================
# MONTE CARLO SIMULATION
# =============================================================================

# Compile Stan model (HEx model variant)
# Note: Requires "model_HEX.stan" file in working directory
mod <- cmdstan_model("model_HEX.stan")

# Run Monte Carlo simulation in parallel
# Using mclapply for parallel processing across CPU cores
# Note: Adjust mc.cores based on available hardware
res <- mclapply(scenarios, function(a) s_replicate(a, mod), mc.cores = 4)

# =============================================================================
# RESULTS ANALYSIS AND VISUALIZATION  
# =============================================================================

#' Precision Analysis for Multiple Intervention Effects
#' 
#' Focuses on the precision (posterior standard deviation) of interaction
#' effects, which are typically the most difficult to estimate precisely
#' and represent the key scientific questions in factorial designs.

# Define interaction effects of interest
# These represent two-way and three-way treatment combinations
x_effects <- c("lOR[4]",    # AB vs None  
               "lOR[5]",    # AC vs None
               "lOR[6]",    # BC vs None
               "lOR[7]")    # ABC vs None

# Extract sample sizes from results
n_ed <- sapply(res, function(x) x[["args"]]["n_ed"])

# Extract posterior standard deviations for interaction effects
res_sd_x <- lapply(res, function(x) {
  x[["model_res"]][variable %in% c(x_effects), .(variable, sd)]
})

# Combine results with sample size information
res_sd_x <- rbindlist(
  Map(function(dt, val) {
    dt[, n_ed := val]
    return(dt)
  }, res_sd_x, n_ed)
)

# Compute average precision across Monte Carlo replications
sum_sd_x <- res_sd_x[, .(sd = mean(sd)), keyby = n_ed]

#' Create Precision vs Sample Size Plot
#' 
#' Visualizes the relationship between number of clusters and precision
#' of interaction effect estimates. The horizontal line at 0.13 represents
#' a precision target (somewhat arbitrary but represents "adequate precision").
#' 

p40 <- ggplot() +
  
  # Reference line for precision target  
  geom_hline(yintercept = 0.13, size = .5, color = "grey80") +
  
  # Individual Monte Carlo results (density and intervals)
  stat_halfeye(data = res_sd_x, 
               aes(group = n_ed, y = sd, x = n_ed, 
                   color = after_stat(as.character(.width))), 
               width = 2, fill = "grey72", 
               point_color = "black", .width = c(.5, .90)) +
  
  # Formatting
  scale_color_manual(values = c("#9a0000", "black"),
                     name = "Interval") +
  scale_x_continuous(limits = c(45, 91), breaks = seq(48, 88, by = 8), name = "number of EDs") +
  scale_y_continuous(limits = c(0.08, 0.20), breaks = c(.13, seq(0.10, .20, by = .05)), 
                     name = "sd of posterior distribution") +
  theme(panel.grid = element_blank(),  
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 12),
        legend.position = "none"
  )

# Display results
print(p40)

# =============================================================================
# INTERPRETATION AND SAMPLE SIZE RECOMMENDATION
# =============================================================================

#' INTERPRETATION OF RESULTS:
#' 
#' The plot shows how posterior precision (standard deviation) decreases
#' as the number of clusters increases. Key considerations:
#' 
#' 1. PRECISION TARGET: The horizontal line at SD = 0.13 represents a
#'    somewhat arbitrary precision target - speficic to our appli aiton. 
#'    You may want to adjust this based on clinically meaningful effect sizes.
#' 
#' 2. DIMINISHING RETURNS: The curve shows diminishing returns - initial
#'    increases in sample size provide large precision gains, but
#'    additional clusters provide smaller improvements.
#' 
#' 3. INTERACTION EFFECTS: This analysis focuses on interaction effects
#'    (lOR[4]-lOR[7]) which are typically hardest to estimate precisely.
#'    Main effects will generally have better precision.
#' 
#' 4. NULL SCENARIO: Results are based on null hypothesis (no true effects).
#'    Precision under alternative hypotheses may differ slightly.
#' 
#' SAMPLE SIZE RECOMMENDATION:
#' Based on the precision target of SD = 0.13, approximately [X] clusters
#' would be needed. However, consider:
#' 
#' - Clinical meaningfulness of the precision target
#' - Cost/feasibility constraints  
#' - Balance between Type I error, power, and precision
#' - Sensitivity to prior specifications
#' 
#' EXTENSIONS:
#' - Test different precision targets
#' - Analyze main effects separately  
#' - Include power analysis for non-null scenarios
#' - Sensitivity analysis across different prior specifications

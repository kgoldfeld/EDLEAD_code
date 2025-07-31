#' ---
#' title: "Factorial Intervention Simulation Study"
#' author: "Keith Goldfeld"
#' date: "`r Sys.Date()`"
#' description: "Monte Carlo simulation study comparing Bayesian and frequentist
#'              approaches for analyzing 2³ factorial cluster randomized trials"
#' ---

# =============================================================================
# REQUIRED LIBRARIES
# =============================================================================

library(cmdstanr)     # Stan interface for Bayesian modeling
library(simstudy)     # Data simulation tools
library(data.table)   # Enhanced data manipulation
library(posterior)    # Posterior analysis tools
library(glmmTMB)      # Generalized linear mixed models

#' Bayesian Model Fitting and Posterior Summarization
#'
#' Fits a Bayesian model to simulated trial data using Stan and extracts
#' posterior summaries for key parameters.
#'
#' @param generated_data A data.table containing simulated trial data
#' @param mod A compiled Stan model object (from cmdstan_model)
#' @param argsvec Named vector containing simulation parameters and hyperparameters
#'
#' @return A data.table with posterior summaries for each parameter:
#'   \itemize{
#'     \item{variable: parameter name (lOR, tau, sigma, delta)}
#'     \item{p.025, p.25, p.50, p.75, p.95, p.975: posterior quantiles}
#'     \item{sd: posterior standard deviation}
#'     \item{var: posterior variance}
#'     \item{p_less_zero: probability parameter is negative}
#'     \item{p_meaningful: probability of meaningful effect (< -0.223, OR < 0.8)}
#'   }
#'
#' @details The function:
#'   \enumerate{
#'     \item Converts data to Stan format using internal dt_to_list helper
#'     \item Fits Bayesian model with robust MCMC settings
#'     \item Extracts posterior draws for treatment effects and model parameters
#'     \item Computes comprehensive posterior summaries
#'   }
#'
#'   MCMC Settings:
#'   \itemize{
#'     \item 4 chains with 500 warmup + 2500 sampling iterations
#'     \item High adapt_delta (0.98) and max_treedepth (20) for complex models
#'     \item Parallel chain execution for efficiency
#'   }
#'
#' @examples
#' \dontrun{
#' posterior_summary <- s_bayes(generated_data, mod, argsvec)
#' }

s_bayes <- function(generated_data, mod, argsvec) {
  
  #' Convert Data.Table to Stan Data List
  #' 
  #' Internal helper function that transforms simulated trial data
  #' into the list format required by Stan models.
  
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                               ## number of observations 
    x_abc <- model.matrix(~r1*r2*r3, data = dx) # Full factorial design matrix
    
    y <- dx[, y]                                # Binary outcomes
    
    N_ED <- dx[, length(unique(ed))]            # Number of clusters
    ed <- dx[, ed]                              # Cluster identifiers
    svals <- c(s1, s2, s3, s4, s5, s6)          # Prior hyperparameters
    
    list(N_ED = N_ED, N = N, x_abc = x_abc, ed = ed, y = y, svals = svals)
  }
  
  # Import simulation parameters into local environment
  list2env(as.list(argsvec), envir = environment())
  
  # Fit Bayesian model with robust MCMC settings
  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,                    # Suppress progress output
    chains = 4L,                    # Multiple chains for convergence
    parallel_chains = 4L,           # Parallel execution
    iter_warmup = 500,              # Warmup iterations
    iter_sampling = 2500,           # Sampling iterations per chain
    adapt_delta = 0.98,             # High acceptance rate for difficult geometry
    max_treedepth = 20,             # Allow deep trees for complex posteriors
    show_messages = FALSE           # Suppress Stan messages
  )
  
  # Convert to posterior draws format
  posterior <- as_draws_array(fit$draws())
  
  # Extract and summarize key parameters
  sumbayes <- as.data.table(posterior)[
    substr(variable, 1, 3) == "lOR" |          # Log-odds ratios (treatment effects)
      substr(variable, 1, 3) == "tau" |        # Regression coefficients
      substr(variable, 1, 5) == "sigma" |      # Scale parameters
      substr(variable, 1, 5) == "delta",       # Additional parameters
    .(
      p.025 = quantile(value, 0.025),          # 95% credible interval bounds
      p.25 = quantile(value, 0.25),            # IQR bounds
      p.50 = quantile(value, 0.50),            # Posterior median
      p.75 = quantile(value, 0.75),            # IQR bounds
      p.95 = quantile(value, 0.95),            # 90% credible interval
      p.975 = quantile(value, 0.975),          # 95% credible interval bounds
      sd = sd(value),                          # Posterior standard deviation
      var = var(value),                        # Posterior variance
      p_less_zero = mean(value < 0),           # Probability of negative effect
      p_meaningful = mean(value < -0.223)      # Probability of meaningful effect
    ), keyby = variable]
  
  return(sumbayes) # model_results is a data.table
  
}

#' Frequentist Analysis Using Generalized Linear Mixed Models
#'
#' Performs hierarchical model selection using likelihood ratio tests to
#' determine the appropriate level of complexity for the factorial design.
#'
#' @param dd A data.table containing the trial data for analysis
#'
#' @return List containing:
#'   \itemize{
#'     \item{glm0: Null model coefficients (intercept only + random effects)}
#'     \item{glm1: Main effects model coefficients}
#'     \item{glm2: Two-way interactions model coefficients}
#'     \item{glm3: Three-way interaction model coefficients}
#'     \item{glmg: Group model coefficients (factor(Grp) + random effects)}
#'     \item{h_conclusion: Hierarchical model selection conclusion}
#'   }
#'
#' @details The function fits five nested models and uses sequential likelihood
#'   ratio tests (LRT) with α = 0.05 to determine model complexity:
#'   \enumerate{
#'     \item Test main effects vs null model
#'     \item If significant, test 2-way interactions vs main effects
#'     \item If significant, test 3-way interaction vs 2-way model
#'   }
#'
#'   Models fitted:
#'   \itemize{
#'     \item{Null: y ~ 1 + (1|ed)}
#'     \item{Main effects: y ~ r1 + r2 + r3 + (1|ed)}
#'     \item{Two-way: y ~ r1*r2 + r1*r3 + r2*r3 + (1|ed)}
#'     \item{Three-way: y ~ r1*r2*r3 + (1|ed)}
#'     \item{Group: y ~ factor(Grp) + (1|ed)}
#'   }
#'
#' @examples
#' \dontrun{
#' freq_results <- s_freq(generated_data)
#' print(freq_results$h_conclusion)
#' }

s_freq <- function(dd) {
  
  # Fit nested models for hierarchical testing
  
  glmfit_3 <- glmmTMB(y ~ r1*r2*r3 + (1|ed), family=binomial, data = dd)
  glmfit_2 <- glmmTMB(y ~ r1*r2 + r1*r3 + r2*r3 + (1|ed), family=binomial, data = dd)
  glmfit_1 <- glmmTMB(y ~ r1 + r2 + r3 + (1|ed), family=binomial, data = dd)
  glmfit_0 <- glmmTMB(y ~ 1 + (1|ed), family=binomial, data = dd)
  glmfit_g <- glmmTMB(y ~ factor(Grp) + (1|ed), family=binomial, data = dd)
  
  # Sequential likelihood ratio tests for model selection
  
  LRT_p1 <- anova(glmfit_1, glmfit_0)$`Pr(>Chisq)`[2]  # Main effects test
  h_conclusion <- NA
  
  if (LRT_p1 <= 0.05) { # main effects present
    
    LRT_p2 <- anova(glmfit_2, glmfit_1)$`Pr(>Chisq)`[2]  # 2-way interactions test
    
    if (LRT_p2 <= 0.05) { # 2-way effects present
      
      LRT_p3 <- anova(glmfit_3, glmfit_2)$`Pr(>Chisq)`[2]  # 3-way interaction test
      if (LRT_p3 <= 0.05) { 
        h_conclusion = "3-way interaction"
      } else {
        h_conclusion = "2-way interaction"
      }
    } else {
      h_conclusion = "main effects only"
    }
  } else {
    h_conclusion = "no effects"
  }
  
  # Return model summaries and hierarchical conclusion
  return(list(
    glm0 = as.data.table(summary(glmfit_0)$coefficients$cond, keep.rownames = "param"),
    glm1 = as.data.table(summary(glmfit_1)$coefficients$cond, keep.rownames = "param"),
    glm2 = as.data.table(summary(glmfit_2)$coefficients$cond, keep.rownames = "param"),
    glm3 = as.data.table(summary(glmfit_3)$coefficients$cond, keep.rownames = "param"),
    glmg = as.data.table(summary(glmfit_g)$coefficients$cond, keep.rownames = "param"),
    h_conclusion = h_conclusion
  ))
  
}

#' Complete Simulation Replication
#'
#' Executes one complete simulation replication including data generation,
#' Bayesian model fitting, and frequentist analysis.
#'
#' @param argsvec Named vector containing all simulation parameters:
#'   \itemize{
#'     \item{n_ed: number of clusters}
#'     \item{icc: intracluster correlation}
#'     \item{n_quarters: number of time periods}
#'     \item{grp1-grp4: cluster size parameters}
#'     \item{t_0, t_a, t_b, t_c: treatment effect parameters}
#'     \item{x_ab, x_ac, x_bc, x_abc: interaction parameters}
#'     \item{s1-s6: prior hyperparameters}
#'   }
#' @param mod List containing two compiled Stan model objects [mod_HEX, mod_nonX]
#'
#' @return List containing:
#'   \itemize{
#'     \item{args: input simulation parameters}
#'     \item{model_Q: posterior summaries from first Bayesian model}
#'     \item{model_S: posterior summaries from second Bayesian model}
#'     \item{model_glm: frequentist analysis results}
#'   }
#'
#' @details This function represents one Monte Carlo iteration:
#'   \enumerate{
#'     \item Generate simulated dataset using s_define() and s_generate()
#'     \item Fit two different Bayesian models
#'     \item Perform frequentist hierarchical model selection
#'     \item Package all results with input parameters
#'   }
#'
#' @examples
#' \dontrun{
#' replication_result <- s_replicate(argsvec, list(mod_HEX, mod_nonX))
#' }

s_replicate <- function(argsvec, mod) {
  
  # Generate simulated dataset
  list_of_defs <- s_define()                        # Define data structure
  generated_data <- s_generate(list_of_defs, argsvec)  # Generate one dataset
  
  # Fit Bayesian models
  model_bayes_1 <- s_bayes(generated_data, mod[[1]], argsvec)  # First model
  model_bayes_2 <- s_bayes(generated_data, mod[[2]], argsvec)  # Second model
  
  # Fit frequentist models
  model_freq <- s_freq(generated_data)
  
  #--- Package summary statistics ---#
  
  summary_stats <- c(
    list(args = argsvec,              # Input parameters
         model_Q =  model_bayes_1,    # First Bayesian model results
         model_S = model_bayes_2,     # Second Bayesian model results
         model_glm = model_freq)      # Frequentist results
  )
  
  return(summary_stats) # summary_stats is a list
}

#' Create Parameter Scenarios for Simulation Study
#'
#' Generates all combinations of specified parameter values for the
#' Monte Carlo simulation study.
#'
#' @param ... Named arguments representing parameter values to combine
#'
#' @return List of parameter vectors, each representing one simulation scenario
#'
#' @details Uses expand.grid() to create all parameter combinations and splits
#'   into list format suitable for parallel processing with mclapply().
#'
#' @examples
#' \dontrun{
#' scenarios <- scenario_list(n_ed = c(48, 64, 80), icc = 0.015)
#' }

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

#--- Study Design Parameters ---#

n_ed <- 8 * 10                        # Number of emergency departments (80)
icc <- 0.015                          # Intracluster correlation coefficient
n_quarters <- 2                       # Number of time periods per cluster

#--- Cluster Size Parameters ---#
# Base values: mean 40 (distribution: .55, .15, .15, .15)

grp1 <- 60                            # Small clusters baseline
grp2 <- 90                            # Medium-small clusters baseline
grp3 <- 190                           # Medium-large clusters baseline
grp4 <- 300                           # Large clusters baseline

### Choose a mean <-------------------

mean <- 40        # possible values = 25, 30, 35, 40
ratio <- mean/40  # Scaling factor from base scenario

# Apply scaling to all cluster size groups
grp1 <- grp1 * ratio
grp2 <- grp2 * ratio
grp3 <- grp3 * ratio
grp4 <- grp4 * ratio

#--- Treatment Effect Parameters ---#
# OR for single intervention: 0.80, OR for 2 interventions: 0.70

t_0 <- -0.40      # Baseline log-odds (approximately 40% control rate)
t_a <- -0.20      # Main effect A (OR ≈ 0.82)
t_b <- 0          # Main effect B (null)
t_c <- -0.20      # Main effect C (OR ≈ 0.82)
x_ab <- 0         # A×B interaction (null)
x_ac <- 0.05      # A×C interaction (small positive)
x_bc <- 0         # B×C interaction (null)
x_abc <- 0.0      # A×B×C three-way interaction (null)

# Alternative null scenario (commented out):
# t_0 <- -0.40
# t_a <- 0
# t_b <- 0
# t_c <- 0
# x_ab <- 0
# x_ac <- 0
# x_bc <- 0
# x_abc <- 0

#--- Prior Hyperparameters ---#

sds <- c(0.3, 0.3, 0.4, 0.4, 0.5, 0.4) # Scenario 1: moderate priors
# sds <- c(0.6, 0.6, 0.8, 0.8, 0.5, 0.8) # Scenario 2: wide priors
# sds <- c(3.0, 0.15, 0.40, 0.40, 3.0, 0.30)  # Scenario 3: inverse gamma

#--- Generate All Parameter Combinations ---#

scenarios <- scenario_list(n_ed = n_ed, icc = icc, n_quarters = n_quarters,
                           grp1 = grp1, grp2 = grp2, grp3 = grp3, grp4 = grp4,
                           t_0 = t_0, t_a = t_a, t_b = t_b, t_c = t_c, 
                           x_ab = x_ab, x_ac = x_ac, x_bc = x_bc, x_abc = x_abc,
                           s1 = sds[1], s2 = sds[2], s3 = sds[3], 
                           s4 = sds[4], s5 = sds[5], s6 = sds[6])

# Replicate each scenario for Monte Carlo simulation
scenarios <- rep(scenarios, each = 1000) # 1000 replications per scenario

# =============================================================================
# MODEL COMPILATION AND SIMULATION EXECUTION
# =============================================================================

#### Compile Stan Models ####

mod_HEX <- cmdstan_model("model_HEX.stan")      # Primary Bayesian model
mod_nonX <- cmdstan_model("model_nonX.stan")    # Alternative Bayesian model
# mod_HEX <- cmdstan_model("model_HEX_ig.stan") # Inverse gamma variant

#### Run Monte Carlo Simulation ####

# Execute simulation in parallel
# Note: Fix typo in original code (mode_nonX should be mod_nonX)
res <- mclapply(scenarios, function(a) s_replicate(a, list(mod_HEX, mod_nonX)), mc.cores = 4)

# =============================================================================
# SIMULATION STUDY NOTES
# =============================================================================

#' STUDY DESIGN NOTES:
#' 
#' This simulation study compares Bayesian and frequentist approaches for
#' analyzing 2³ factorial cluster randomized trials with the following features:
#' 
#' 1. DESIGN: 2³ factorial (8 treatment combinations) with 80 clusters
#' 2. OUTCOME: Binary outcome with logistic regression
#' 3. CLUSTERING: Random intercepts for emergency departments
#' 4. BAYESIAN MODELS: Two different Stan implementations
#' 5. FREQUENTIST: Hierarchical model selection via LRT
#' 
#' KEY RESEARCH QUESTIONS:
#' - Type I error rates and power comparison
#' - Robustness to prior specifications
#' - Model selection performance
#' - Estimation accuracy for interaction effects
#' 
#' COMPUTATIONAL REQUIREMENTS:
#' - 1000 replications × 2 Bayesian models = 2000 MCMC runs
#' - Substantial memory and time requirements
#' - Parallel processing recommended (mc.cores = 4)
#' 
#' DEPENDENCIES:
#' - Requires s_define() and s_generate() functions (not included)
#' - Requires model_HEX.stan and model_nonX.stan files
#' - Original code has typo: "mode_nonX" should be "mod_nonX"

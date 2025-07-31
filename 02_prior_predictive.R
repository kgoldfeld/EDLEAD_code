#' ---
#' title: "Bayesian Prior Predictive Analysis for Cluster Randomized Trial"
#' author: "Keith Goldfeld"
#' ---

# =============================================================================
# REQUIRED LIBRARIES
# =============================================================================

library(cmdstanr)     # Bayesian modeling with Stan
library(ggplot2)      # Data visualization
library(ggpubr)       # Publication-ready plots and plot arrangements
library(data.table)   # Efficient data manipulation (from previous code)

#' Convert Data Table to Stan Data List
#' 
#' Transforms a data.table containing trial data into a list format
#' suitable for Stan model fitting.
#' 
#' @param dx A data.table containing the simulated trial data with columns:
#'   \itemize{
#'     \item{ed: cluster/site identifier}
#'     \item{r1, r2, r3: binary treatment indicators}
#'     \item{y: binary outcome variable}
#'   }
#' @param argsvec Named vector containing prior hyperparameter values:
#'   \itemize{
#'     \item{sigma_tau_m: prior SD for main effect coefficients}
#'     \item{sigma_tau_x: prior SD for interaction coefficients}
#'     \item{ss_m, ss_x, ss_3: hierarchical prior parameters}
#'     \item{ss_ed: prior SD for cluster random effects}
#'   }
#' 
#' @return A list containing Stan data components:
#'   \itemize{
#'     \item{N: number of observations}
#'     \item{N_ED: number of clusters}
#'     \item{x_abc: design matrix with main effects and interactions}
#'     \item{ed: cluster identifiers}
#'     \item{y: binary outcomes}
#'     \item{sigma_*: prior hyperparameters}
#'   }
#' 
#' @details Creates a full factorial design matrix including:
#'   - Intercept
#'   - Main effects (r1, r2, r3)
#'   - Two-way interactions (r1:r2, r1:r3, r2:r3)
#'   - Three-way interaction (r1:r2:r3)

dt_to_list <- function(dx, rep_scenario) {
  
  # Import hyperparameters into local environment
  list2env(as.list(rep_scenario), envir = environment())

  # Extract data dimensions and components
  N <- nrow(dx)                               # Total number of observations 
  x_abc <- model.matrix(~r1*r2*r3, data = dx) # Full factorial design matrix
  y <- dx[, y]                                # Binary outcomes
  N_ED <- dx[, length(unique(ed))]            # Number of clusters
  ed <- dx[, ed]                              # Cluster identifiers
  
  # Return Stan data list with data and hyperparameters
  list(N_ED = N_ED, 
       N = N, 
       x_abc = x_abc, 
       ed = ed, 
       y = y,
       sigma_tau_m = sigma_tau_m,    # Prior SD for main effects
       sigma_tau_x = sigma_tau_x,    # Prior SD for interactions
       sigma_sigma_m = ss_m,         # Hierarchical prior parameter
       sigma_sigma_x = ss_x,         # Hierarchical prior parameter
       sigma_3 = ss_3,               # Hierarchical prior parameter
       sigma_sigma_ed = ss_ed)       # Prior SD for cluster effects
}

#' Generate Prior Predictive Check Plot
#' 
#' Runs prior predictive simulation and creates histogram plots
#' showing the distribution of predicted probabilities for each
#' treatment arm combination.
#' 
#' @param argsvec Named vector containing hyperparameter values for the prior
#'   distributions. Must include sigma_tau_m for panel labeling.
#' 
#' @return A ggplot object showing histograms of predicted probabilities
#'   faceted by treatment arm, with panel title indicating hyperparameter values.
#' 
#' @details The function:
#'   1. Samples from prior distributions using Stan
#'   2. Generates predicted outcomes (y_rep) for each observation
#'   3. Calculates mean predicted probability for each treatment arm
#'   4. Creates histograms showing distribution across prior draws
#'   
#'   Treatment arms are labeled as:
#'   - None: no interventions
#'   - A, B, C: single interventions
#'   - AB, AC, BC: two-way combinations
#'   - ABC: all three interventions
#' 
#' @note Requires a compiled Stan model object 'mod_prior' and 
#'   simulated data 'generated_data' to be available in the global environment.

replicate <- function(rep_scenario) {
  
  # Import hyperparameters into local environment
  list2env(as.list(rep_scenario), envir = environment())
  
  # Create panel labels based on main effect prior SD
  if (sigma_tau_m == 0.4) P <- "Panel 1: "
  if (sigma_tau_m == 0.8) P <- "Panel 2: "
  if (sigma_tau_m == 2.5) P <- "Panel 3: "
  
  # Run prior predictive sampling
  fit_prior <- mod_prior$sample(
    data = dt_to_list(generated_data, rep_scenario),
    chains = 1,                    # Single chain for prior predictive
    iter_sampling = 10000,         # Number of prior draws
    fixed_param = TRUE,            # Prior predictive mode
    refresh = 1000                 # Progress update frequency
  )
  
  # Extract posterior predictive draws
  y_rep_draws <- fit_prior$draws(c("y_rep"))
  
  # Calculate mean predicted probability for each treatment arm
  probs <- rbindlist(lapply(0:7, function(i) {
    vars_of_interest <- which(generated_data$arm == i)
    arm_probs <- apply(y_rep_draws[,,vars_of_interest], 1, mean)
    data.table(arm_probs)
  }), idcol = TRUE)
  
  # Create treatment arm labels
  labels <- c("None","A", "B", "C", "AB", "AC", "BC", "ABC")
  probs[, arm_label := factor(.id, levels = 1:8, labels = labels)]
  
  # Reverse factor levels for display (None at top, ABC at bottom)
  probs[, arm_label := factor(arm_label, levels = rev(labels))]
  
  # Create histogram plot
  ggplot(data = probs, aes(x = arm_probs)) +
    geom_histogram(boundary = 0, 
                   color = "black", 
                   fill = "tomato",  
                   aes(y = after_stat(density)), 
                   bins = 30) +
    facet_grid(arm_label ~ .) +
    ggtitle(bquote(bold(.(P)) ~ nu[tau] == .(sigma_tau_m) * "," ~ 
                     nu[sigma] == .(ss_m))) +
    labs(x = "predicted probability",
         y = "density") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank())
}

# =============================================================================
# COMPILE STAN MODEL
# =============================================================================

# Compile the prior-only model
# Note: Requires "prior_predictive_model.stan" file in working directory
mod_prior <- cmdstan_model("prior_predictive_model.stan")

# =============================================================================
# GENERATE SIMULATED DATASET
# =============================================================================

# Study design parameters
n_ed <- 48                        # Number of clusters (sites)
icc <- 0.015                      # Intracluster correlation coefficient
n_quarters <- 2                   # Number of time periods per cluster
n_quarter <- c(60, 90, 190, 300)  # Expected observations per quarter by group

# Model parameters (under null hypothesis)
t_0 <- -0.85                 # Intercept (log-odds scale)
t_a <- t_b <- t_c <- 0       # Main effect parameters
x_ab <- x_ac <- x_bc <- 0    # Two-way interaction parameters
x_abc <- 0                   # Three-way interaction parameter

# Generate parameter combinations for simulation study
# Note: Requires scenario_list() function from previous code
argsvec <- scenario_list(
  n_ed = n_ed,
  icc = icc,
  n_quarters = n_quarters,
  grp1 = n_quarter[1],       # Small clusters
  grp2 = n_quarter[2],       # Medium-small clusters
  grp3 = n_quarter[3],       # Medium-large clusters
  grp4 = n_quarter[4],       # Large clusters
  t_0 = t_0,
  t_a = t_a,
  t_b = t_b,
  t_c = t_c,
  x_ab = x_ab,
  x_ac = x_ac,
  x_bc = x_bc,
  x_abc = x_abc
)

# Create data definition structure
# Note: Requires s_define() function from previous code
list_of_defs <- s_define()

# Generate single simulated dataset
# Note: Requires s_generate() function from previous code
generated_data <- lapply(argsvec, function(args) {
  s_generate(list_of_defs, args)
})[[1]]

# =============================================================================
# PRIOR SENSITIVITY ANALYSIS SCENARIOS
# =============================================================================

#' Prior Hyperparameter Scenarios
#' 
#' Three scenarios representing different levels of prior informativeness:
#' 
#' Scenario 1 (Informative): Small prior SDs suggesting modest treatment effects
#' Scenario 2 (Moderate): Medium prior SDs allowing moderate treatment effects  
#' Scenario 3 (Weakly informative): Large prior SDs allowing substantial effects
#' 
#' Each scenario specifies:
#' - sigma_tau_m/x: Prior SDs for main effects and interactions
#' - ss_m/x: Hierarchical prior parameters for effect heterogeneity
#' - ss_3: Additional hierarchical parameter
#' - ss_ed: Prior SD for cluster random effects

scenarios <- list(
  # Scenario 1: Informative priors (small treatment effects expected)
  c(sigma_tau_m = 0.4, sigma_tau_x = 0.4, ss_m = 0.3, ss_x = 0.3, 
    ss_3 = 0.3, ss_ed = 0.5),
  
  # Scenario 2: Moderately informative priors  
  c(sigma_tau_m = 0.8, sigma_tau_x = 0.8, ss_m = 0.6, ss_x = 0.6, 
    ss_3 = 0.3, ss_ed = 0.5),
  
  # Scenario 3: Weakly informative priors (large effects possible)
  c(sigma_tau_m = 2.5, sigma_tau_x = 2.5, ss_m = 2.5, ss_x = 2.5, 
    ss_3 = 0.3, ss_ed = 0.5)
)

# =============================================================================
# GENERATE PRIOR PREDICTIVE CHECK PLOTS
# =============================================================================

# Generate plots for each prior scenario
plots <- lapply(scenarios, function(x) replicate(x))

# Combine plots into single figure
p1 <- ggarrange(plotlist = plots, nrow = 1)

# Display the combined plot
print(p1)

# =============================================================================
# INTERPRETATION NOTES
# =============================================================================

# The resulting plots show:
# 
# 1. How different prior specifications affect predicted outcome probabilities
# 2. Whether priors lead to reasonable predictions across treatment arms
# 3. The range of plausible effect sizes implied by each prior
# 
# Key considerations:
# - Do the predicted probabilities span a reasonable range?
# - Are there concerning patterns (e.g., all arms having identical distributions)?
# - Do the priors allow for clinically meaningful effect sizes?
# 
# This analysis helps validate prior choices before fitting to real data.

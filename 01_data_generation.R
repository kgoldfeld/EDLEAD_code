#' ---
#' title: "Generating Data for Simulation Study of Cluster Randomized Trial"
#' author: "Keith Goldfeld"
#' ---

# Load required libraries

library(simstudy)  # For data simulation
library(data.table) # For efficient data manipulation

#' Define Data Generation Structure
#' 
#' Creates the data definition objects for the simulation study.
#' This function sets up:
#' - Cluster-level random effects and group assignments
#' - Conditional quarter definitions based on group membership
#' - Period-level observation counts
#' - Binary outcome model with main effects and interactions
#' 
#' @return A list containing four data definition objects:
#'   \item{d1}{Cluster-level definitions (random effects, group assignment)}
#'   \item{d2}{Period-level definitions (number of observations per period)}
#'   \item{def_qn}{Conditional definitions for quarter sizes by group}
#'   \item{defY}{Binary outcome model with logistic regression structure}
#'   
#' @details The outcome model includes:
#'   - Intercept (t_0)
#'   - Cluster random effect (a)
#'   - Main effects for three factors (t_a, t_b, t_c)
#'   - Two-way interactions (x_ab, x_ac, x_bc)
#'   - Three-way interaction (x_abc)

s_define <- function() {
  
  # Cluster-level data definitions
  # - Random effect 'a' with variance determined by ICC
  # - Group assignment with specified probabilities (55%, 15%, 15%, 15%)
  d1 <- defData(varname = "a", 
                formula = 0, 
                variance = "..revar", 
                dist = "normal", 
                id = "ed") |>
    defData(varname = "s_grp", 
            formula = "0.55;0.15;0.15;0.15", 
            dist = "categorical")
  
  # Conditional definitions for quarter sizes by group
  # Each group has different expected number of observations per quarter
  def_qn <- defCondition(condition = "s_grp == 1", 
                         formula = "..grp1", 
                         dist = "poisson") |>
    defCondition(condition = "s_grp == 2", 
                 formula = "..grp2", 
                 dist = "poisson") |>
    defCondition(condition = "s_grp == 3", 
                 formula = "..grp3", 
                 dist = "poisson") |>
    defCondition(condition = "s_grp == 4", 
                 formula = "..grp4", 
                 dist = "poisson")
  
  # Period-level data: number of observations per period
  d2 <- defDataAdd(varname = "nper", 
                   formula = "quarter_n", 
                   dist = "poisson")
  
  # Binary outcome model with logistic link
  # Includes main effects, two-way interactions, and three-way interaction
  defY <- 
    defDataAdd(varname = "arm", 
      formula = "r1 + r2*2 + r3*4", 
      dist = "nonrandom") |>
    defDataAdd(varname = "y", 
      formula = "..t_0 + a + ..t_a*r1 + ..t_b*r2 + ..t_c*r3 + 
      ..x_ab*r1*r2 + ..x_ac*r1*r3 + ..x_bc*r2*r3 + 
      ..x_abc*r1*r2*r3", 
      dist = "binary", 
      link = "logit")
  
  return(list(d1 = d1, 
              d2 = d2, 
              def_qn = def_qn, 
              defY = defY))
}

#' Generate Simulated Dataset
#' 
#' Generates a single simulated dataset based on the provided definitions
#' and parameter values.
#' 
#' @param list_of_defs List of data definition objects from s_define()
#' @param argsvec Named vector of parameter values for the simulation
#' 
#' @return A data.table containing the simulated dataset with:
#'   - Cluster identifiers (ed)
#'   - Time period identifiers
#'   - Treatment assignments (r1, r2, r3)
#'   - Individual observations (id)
#'   - Binary outcomes (y)
#'   
#' @details The function:
#'   1. Generates clusters with random effects and group assignments
#'   2. Assigns quarter sizes based on group membership
#'   3. Performs stratified randomization by cluster size
#'   4. Creates longitudinal structure with multiple periods
#'   5. Generates individual-level observations within periods
#'   6. Simulates binary outcomes based on the specified model

s_generate <- function(list_of_defs, argsvec) {
  
  # Import definitions and parameters into local environment
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  # Calculate random effect variance based on ICC
  revar <- iccRE(icc, dist = "binary")
  
  # Generate cluster-level data
  dd <- genData(n_ed, d1)
  
  # Add conditional quarter definitions
  dd <- addCondition(def_qn, dd, "quarter_n")
  
  # Stratified randomization by cluster size group
  # Ensures balanced treatment assignment within each size stratum
  dd <- rbindlist(
    lapply(c(1:4), function(x) {
      addMultiFac(dd[s_grp == x], 
                  nFactors = 3, 
                  colNames = c("r1", "r2", "r3"))
    })
  )
  setkey(dd, "ed")
  
  # Add longitudinal structure (multiple quarters per cluster)
  dd <- addPeriods(dd, nPeriods = n_quarters, idvars = "ed")
  
  # Add number of observations per period
  dd <- addColumns(d2, dd)
  
  # Generate individual-level observations within each period
  dd <- genCluster(dd, "timeID", "nper", "id")
  
  # Generate binary outcomes
  dd <- addColumns(defY, dd)
  
  # Alternative grouping
  
  dd[, Grp := arm + 1]
  dd[arm == 3, Grp := 5]
  dd[arm == 4, Grp := 4]
  dd[, Grp := factor(Grp)]
  
  return(dd)
}

#' Create Parameter Scenarios
#' 
#' Generates all combinations of specified parameter values for simulation study.
#' 
#' @param ... Named arguments specifying parameter values to vary
#' 
#' @return List of parameter vectors, each representing one simulation scenario
#' 
#' @examples
#' scenarios <- scenario_list(n_clusters = c(20, 40), icc = c(0.01, 0.05))

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Study design parameters
n_ed <- c(40, 48, 56)             # Number of clusters (sites)
icc <- 0.015                      # Intracluster correlation coefficient
n_quarters <- 2                   # Number of time periods per cluster
n_quarter <- c(60, 90, 190, 300)  # Expected observations per quarter by group

# Model parameters
t_0 <- -0.40                 # Intercept (log-odds scale)
t_a <- t_b <- t_c <- 0       # Main effect parameters (null hypothesis)
x_ab <- x_ac <- x_bc <- 0    # Two-way interaction parameters (null)
x_abc <- 0                   # Three-way interaction parameter (null)

# Generate all parameter combinations for simulation study
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

# =============================================================================
# GENERATE SIMULATED DATASETS
# =============================================================================

# Create data definition structure
list_of_defs <- s_define()

# Generate datasets for all parameter scenarios
# Each element of the list contains one simulated dataset
generated_data <- lapply(argsvec, function(args) {
  s_generate(list_of_defs, args)
})


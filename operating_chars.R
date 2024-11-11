library(cmdstanr)
library(simstudy)
library(data.table)
library(posterior)

## Need to specify mean, model hierarchy, and output name

s_define <- function() {
  
  #--- data definition code ---#
  
  d1 <- defData(varname = "a", formula = 0, variance = "..revar", dist = "normal", id="ed")
  d1 <- defData(d1, varname = "s_grp", 
    formula = "0.55;0.15;0.15;0.15", dist = "categorical")
  
  d2 <- defDataAdd(varname = "nper", formula = "quarter_n", dist="poisson")
  
  defS <- defCondition(condition = "s_grp == 1",  formula = "..grp1", dist = "poisson")
  defS <- defCondition(defS, condition = "s_grp == 2",  formula = "..grp2", dist = "poisson")
  defS <- defCondition(defS, condition = "s_grp == 3",  formula = "..grp3", dist = "poisson")
  defS <- defCondition(defS, condition = "s_grp == 4",  formula = "..grp4", dist = "poisson")
  
  defY <- defDataAdd(varname = "y", 
    formula = "..t_0 + a + ..t_a*r1 + ..t_b*r2 + ..t_c*r3 + ..x_ab*r1*r2 + ..x_ac*r1*r3 + ..x_bc*r2*r3 + ..x_abc*r1*r2*r3",
    dist = "binary", link="logit")
  
  return(list(d1 = d1, d2 = d2, defS = defS, defY = defY)) 
  
}

s_generate <- function(list_of_defs, argsvec) {
  
  list2env(list_of_defs, envir = environment())
  list2env(as.list(argsvec), envir = environment())
  
  revar <- iccRE(icc, dist = "binary")
  print(revar)
  
  dd <- genData(n_ed, d1)
  dd <- addCondition(defS, dd, "quarter_n")
  
  dd <- rbindlist(lapply(c(1:4), 
    function(x) addMultiFac(dd[s_grp == x], nFactors = 3, colNames = c("r1", "r2", "r3")))
  )
  setkey(dd, "ed")
  
  dd <- addPeriods(dd, nPeriods = n_quarters, idvars = "ed")
  dd <- addColumns(d2, dd)
  dd <- genCluster(dd, "timeID", "nper", "id")
  dd <- addColumns(defY, dd)
  
  print("done")
  
  return(dd)
  
}

s_model <- function(generated_data, mod) {
  
  dt_to_list <- function(dx) {
    
    N <- nrow(dx)                               ## number of observations 
    x_abc <- model.matrix(~r1*r2*r3, data = dx)
    
    y <- dx[, y]
    
    N_ED <- dx[, length(unique(ed))]
    ed <- dx[, ed]
    
    list(N_ED = N_ED, N = N, x_abc = x_abc, ed = ed, y = y)
  }
  
  print("starting fit")
  
  fit <- mod$sample(
    data = dt_to_list(generated_data),
    refresh = 0,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = 500,
    iter_sampling = 2500,
    adapt_delta = 0.98,
    max_treedepth = 20,
    show_messages = FALSE
  )
  
  print("ending fit")
  
  posterior <- as_draws_array(fit$draws())
  
  print("preparing final results")
  
  sumstats <- as.data.table(posterior)[
    substr(variable, 1, 3) == "lOR" |
    substr(variable, 1, 5) == "sigma" |
    substr(variable, 1, 5) == "delta", 
    .(
      p.025 = quantile(value, 0.025),
      p.25 = quantile(value, 0.25),
      p.50 = quantile(value, 0.50),
      p.75 = quantile(value, 0.75),
      p.95 = quantile(value, 0.95),
      p.975 = quantile(value, 0.975),
      sd = sd(value),
      var = var(value),
      p_less_zero = mean(value < 0),
      p_meaningful = mean(value < -0.223)
    ), keyby = variable]
  
  return(sumstats) # model_results is a data.table
  
}

s_replicate <- function(argsvec, mod) {
  
  set_cmdstan_path(path = "/gpfs/share/apps/cmdstan/2.34.1")
  
  list_of_defs <- s_define()
  generated_data <- s_generate(list_of_defs, argsvec)
  model_results_1 <- s_model(generated_data, mod[[1]])
  model_results_2 <- s_model(generated_data, mod[[2]])
  
  #--- summary statistics ---#
  
  summary_stats <- list(
    model1 = data.table(t(argsvec), model_results_1),
    model2 = data.table(t(argsvec), model_results_2)
  )
  
  return(summary_stats) # summary_stats is a data.table
}

#--- Set arguments ---#

scenario_list <- function(...) {
  argmat <- expand.grid(...)
  return(asplit(argmat, MARGIN = 1))
}

n_ed <- 8 * c(10)
icc <- c(0.005, 0.010, 0.015, 0.020)

n_quarters <- 2

# base values: mean 40 (distribution: .55, .15, .15, .15)

grp1 <- 60
grp2 <- 90
grp3 <- 190
grp4 <- 300

### Choose a mean <-------------------

mean <- c(25, 30, 35, 40, 45)        # possible values = 25, 30, 35, 40
ratio <- mean/40

grp1 <- grp1 * ratio
grp2 <- grp2 * ratio
grp3 <- grp3 * ratio
grp4 <- grp4 * ratio

# OR for single intervention: 0.80, OR for 2 interventions: 0.70

# t_0 <- -0.85
# t_a <- -0.20
# t_b <- 0
# t_c <- -0.10
# x_ab <- 0
# x_ac <- 0.05
# x_bc <- 0
# x_abc <- 0.0

t_0 <- -.85
t_a <- 0
t_b <- 0
t_c <- 0
x_ab <- 0
x_ac <- 0
x_bc <- 0
x_abc <- 0

scenarios <- lapply(1:4, 
  function(i) {
    a <- scenario_list(n_ed = n_ed, icc = icc[i], n_quarters = n_quarters,
          grp1 = grp1, grp2 = grp2, grp3 = grp3, grp4 = grp4,
          t_0 = t_0, t_a = t_a, t_b = t_b, t_c = t_c, 
          x_ab = x_ab, x_ac = x_ac, x_bc = x_bc, x_abc = x_abc)
    a[[1]]}
)

scenarios <- rep(scenarios, each = 400)

mod1 <- cmdstan_model("bayes_hpc_1level.stan")
mod2 <- cmdstan_model("bayes_hpc_2level.stan")

summary_stats <- lapply(scenarios, function(a) s_replicate(a, list(mod1, mod2)))




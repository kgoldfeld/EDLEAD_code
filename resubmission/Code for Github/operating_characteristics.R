  library(cmdstanr)
  library(simstudy)
  library(data.table)
  library(posterior)
  library(glmmTMB)
  
  s_bayes <- function(generated_data, mod, argsvec) {
    
    dt_to_list <- function(dx) {
      
      N <- nrow(dx)                               ## number of observations 
      x_abc <- model.matrix(~r1*r2*r3, data = dx)
      
      y <- dx[, y]
      
      N_ED <- dx[, length(unique(ed))]
      ed <- dx[, ed]
      svals <- c(s1, s2, s3, s4, s5, s6)
      
      list(N_ED = N_ED, N = N, x_abc = x_abc, ed = ed, y = y, svals = svals)
    }
    
    list2env(as.list(argsvec), envir = environment())
    
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
    
    posterior <- as_draws_array(fit$draws())
    
    sumbayes <- as.data.table(posterior)[
      substr(variable, 1, 3) == "lOR" |
        substr(variable, 1, 3) == "tau" |
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
    
    return(sumbayes) # model_results is a data.table
    
  }
  
  s_freq <- function(dd) {
    
    # Fit models
    
    glmfit_3 <- glmmTMB(y ~ r1*r2*r3 + (1|ed), family=binomial, data = dd)
    glmfit_2 <- glmmTMB(y ~ r1*r2 + r1*r3 + r2*r3 + (1|ed), family=binomial, data = dd)
    glmfit_1 <- glmmTMB(y ~ r1 + r2 + r3 + (1|ed), family=binomial, data = dd)
    glmfit_0 <- glmmTMB(y ~ 1 + (1|ed), family=binomial, data = dd)
    glmfit_g <- glmmTMB(y ~ factor(Grp) + (1|ed), family=binomial, data = dd)
    
    # Global test for main effects
    
    LRT_p1 <- anova(glmfit_1, glmfit_0)$`Pr(>Chisq)`[2]
    h_conclusion <- NA
    
    if (LRT_p1 <= 0.05) { # main effects present
      
      LRT_p2 <- anova(glmfit_2, glmfit_1)$`Pr(>Chisq)`[2]
      
      if (LRT_p2 <= 0.05) { # 2-way effects present
        
        LRT_p3 <- anova(glmfit_3, glmfit_2)$`Pr(>Chisq)`[2]
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
    
    return(list(
      glm0 = as.data.table(summary(glmfit_0)$coefficients$cond, keep.rownames = "param"),
      glm1 = as.data.table(summary(glmfit_1)$coefficients$cond, keep.rownames = "param"),
      glm2 = as.data.table(summary(glmfit_2)$coefficients$cond, keep.rownames = "param"),
      glm3 = as.data.table(summary(glmfit_3)$coefficients$cond, keep.rownames = "param"),
      glmg = as.data.table(summary(glmfit_g)$coefficients$cond, keep.rownames = "param"),
      h_conclusion = h_conclusion
    ))
    
  }
  
  s_replicate <- function(argsvec, mod) {
    
    list_of_defs <- s_define()
    generated_data <- s_generate(list_of_defs, argsvec)
    model_bayes_1 <- s_bayes(generated_data, mod[[1]], argsvec)
    model_bayes_2 <- s_bayes(generated_data, mod[[2]], argsvec)
    model_freq <- s_freq(generated_data)
    
    #--- summary statistics ---#
    
    summary_stats <- c(
      list(args = argsvec, model_Q =  model_bayes_1, model_S = model_bayes_2,
      model_glm = model_freq)
    )
    
    return(summary_stats) # summary_stats is a data.table
  }
  
  #--- Set arguments ---#
  
  scenario_list <- function(...) {
    argmat <- expand.grid(...)
    return(asplit(argmat, MARGIN = 1))
  }
  
  
  n_ed <- 8 * 10
  icc <- 0.015
  
  n_quarters <- 2
  
  # base values: mean 40 (distribution: .55, .15, .15, .15)
  
  grp1 <- 60
  grp2 <- 90
  grp3 <- 190
  grp4 <- 300
  
  ### Choose a mean <-------------------
  
  mean <- 40        # possible values = 25, 30, 35, 40
  ratio <- mean/40
  
  grp1 <- grp1 * ratio
  grp2 <- grp2 * ratio
  grp3 <- grp3 * ratio
  grp4 <- grp4 * ratio
  
  # OR for single intervention: 0.80, OR for 2 interventions: 0.70
  
  t_0 <- -0.40
  t_a <- -0.20
  t_b <- 0
  t_c <- -0.20
  x_ab <- 0
  x_ac <- 0.05
  x_bc <- 0
  x_abc <- 0.0
  
  # t_0 <- -0.40
  # t_a <- 0
  # t_b <- 0
  # t_c <- 0
  # x_ab <- 0
  # x_ac <- 0
  # x_bc <- 0
  # x_abc <- 0
  
  sds <- c(0.3, 0.3, 0.4, 0.4, 0.5, 0.4) # 1
  # sds <- c(0.6, 0.6, 0.8, 0.8, 0.5, 0.8) # 2
  # sds <- c(3.0, 0.15, 0.40, 0.40, 3.0, 0.30)  # 3: inverse gamma
  
  scenarios <- scenario_list(n_ed = n_ed, icc = icc, n_quarters = n_quarters,
                grp1 = grp1, grp2 = grp2, grp3 = grp3, grp4 = grp4,
                t_0 = t_0, t_a = t_a, t_b = t_b, t_c = t_c, 
                x_ab = x_ab, x_ac = x_ac, x_bc = x_bc, x_abc = x_abc,
                s1 = sds[1], s2 = sds[2], s3 = sds[3], s4 = sds[4], s5 = sds[5], s6 = sds[6])
  
  scenarios <- rep(scenarios, each = 1000) # for final change back to 1000
  
  ####
  
  print("compling stan")
  
  mod_HEX <- cmdstan_model("model_HEX.stan")
  mod_nonX <- cmdstan_model("model_nonX.stan")
  # mod_HEX <- cmdstan_model("model_HEX_ig.stan")
  
  res <- mclapply(scenarios, function(a) s_replicate(a, list(mod_HEX, mode_nonX)), mc.cores = 4)
  





# Bayesian Analysis of Factorial Cluster Randomized Trials

A comprehensive simulation study comparing Bayesian and frequentist approaches for analyzing 2³ factorial cluster randomized trials, with applications to emergency department interventions.

## Overview

This repository contains code for conducting Monte Carlo simulation studies to evaluate statistical methods for factorial cluster randomized trials. The framework is designed for healthcare settings where multiple interventions are tested simultaneously across clustered units (e.g., emergency departments).

## Key Features

- **2³ Factorial Design**: Analyzes trials with three binary interventions (A, B, C) and all interactions
- **Cluster Randomization**: Accounts for clustering within healthcare sites using random effects
- **Bayesian Methods**: Implements hierarchical models with different prior specifications
- **Frequentist Comparison**: Includes hierarchical model selection via likelihood ratio tests
- **Power Analysis**: Determines required sample sizes based on posterior precision criteria
- **Prior Sensitivity**: Evaluates robustness to different prior assumptions

## Repository Structure

### R Programs

| File | Description |
|------|-------------|
| `01_data_generation.R` | Core functions for simulating factorial cluster randomized trial data |
| `02_prior_predictive.R` | Prior predictive checks to validate prior specifications |
| `03_bayesian_power.R` | Monte Carlo simulation for Bayesian power analysis and sample size determination |
| `04_comparison_study.R` | Comprehensive comparison of Bayesian and frequentist approaches |

### Stan Models

| File | Description |
|------|-------------|
| `model_HEX.stan` | **Hierarchical Exchangeable (HEx) Model**: Main Bayesian model with hierarchical priors for treatment effects |
| `model_nonX.stan` | **Non-Exchangeable Model**: Alternative specification without hierarchical structure |
| `model_HEX_ig.stan` | **HEx with Inverse Gamma Priors**: Variant using inverse gamma priors for scale parameters |
| `prior_predictive_model.stan` | **Prior Predictive Model**: Generates samples from prior distributions for prior checking |

## Study Design

### Treatment Structure
- **Control**: No interventions
- **Single Interventions**: A only, B only, C only
- **Two-way Combinations**: AB, AC, BC
- **Three-way Combination**: ABC

### Data Generation Process
1. **Cluster Level**: Random effects and group assignments (4 size categories)
2. **Time Structure**: Multiple quarters per cluster
3. **Individual Level**: Binary outcomes with logistic regression
4. **Effect Structure**: Main effects, two-way interactions, three-way interaction

### Model Specifications

The Bayesian models use the following hierarchical structure:

```
logit(P(y_i = 1)) = α_j[i] + τ_a*A_i + τ_b*B_i + τ_c*C_i + 
                    τ_ab*A_i*B_i + τ_ac*A_i*C_i + τ_bc*B_i*C_i + τ_abc*A_i*B_i*C_i
```

Where:
- `α_j[i]`: Random effect for cluster j
- `τ_a, τ_b, τ_c`: Main effects for treatments A, B, C
- `τ_ab, τ_ac, τ_bc`: Two-way interaction effects  
- `τ_abc`: Three-way interaction effect

## Getting Started

### Prerequisites

```r
# Required R packages
library(simstudy)     # Data simulation
library(data.table)   # Data manipulation
library(cmdstanr)     # Stan interface
library(posterior)    # Posterior analysis
library(glmmTMB)      # Frequentist mixed models
library(ggplot2)      # Visualization
library(ggpubr)       # Plot arrangements
library(parallel)     # Parallel computing
```

### Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/factorial-cluster-trials.git
cd factorial-cluster-trials
```

2. Install required R packages:
```r
install.packages(c("simstudy", "data.table", "glmmTMB", "ggplot2", "ggpubr", "parallel"))
```

3. Install cmdstanr and Stan:
```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

### Basic Usage

#### 1. Data Generation
```r
source("01_data_generation.R")

# Define simulation parameters
list_of_defs <- s_define()
argsvec <- c(n_ed = 48, icc = 0.015, n_quarters = 2, ...)

# Generate simulated dataset
generated_data <- s_generate(list_of_defs, argsvec)
```

#### 2. Prior Predictive Check
```r
source("02_prior_predictive.R")

# Compile prior predictive model
mod_prior <- cmdstan_model("prior_predictive_model.stan")

# Generate prior predictive plots
plots <- lapply(scenarios, function(x) replicate(x))
```

#### 3. Bayesian Analysis
```r
source("03_bayesian_power.R")

# Compile main model
mod <- cmdstan_model("model_HEX.stan")

# Run power analysis
results <- mclapply(scenarios, function(a) s_replicate(a, mod), mc.cores = 4)
```

## Simulation Parameters

### Study Design Parameters
- **Number of clusters**: 48-88 emergency departments
- **Intracluster correlation**: 0.015
- **Time periods**: 2 quarters per cluster
- **Cluster sizes**: Variable (60-300 patients per quarter)

### Treatment Effect Parameters
- **Baseline log-odds**: -0.40 (≈40% control rate)
- **Single intervention OR**: 0.80-0.82
- **Interaction effects**: Small positive or null

### Prior Specifications
Three scenarios representing different levels of informativeness:
1. **Informative**: Small prior SDs (σ = 0.3-0.4)
2. **Moderate**: Medium prior SDs (σ = 0.6-0.8)  
3. **Weakly informative**: Large prior SDs (σ = 2.5)

## Key Analyses

### 1. Prior Predictive Checking
- Validates prior assumptions by examining implied outcome distributions
- Ensures priors allow for clinically meaningful effect sizes
- Identifies potential prior-data conflicts

### 2. Power Analysis
- Determines required sample size for adequate posterior precision
- Target: Posterior SD ≤ 0.13 for interaction effects
- Monte Carlo estimation with 1000 replications per scenario

### 3. Method Comparison
- Bayesian hierarchical models vs. frequentist hierarchical model selection
- Type I error rates and power comparison
- Robustness to prior specifications

## Computational Requirements

- **Memory**: 8+ GB RAM recommended
- **CPU**: Multi-core processor for parallel processing
- **Time**: Full simulation study requires several hours/days on HPC
- **Dependencies**: Stan installation required

For large-scale simulations, high-performance computing (HPC) is recommended. Contact the authors for HPC setup guidance.

## Key Results

### Sample Size Recommendations
Based on posterior precision targets:
- **Interaction effects**: ~72-80 clusters needed for SD ≤ 0.13
- **Main effects**: Smaller sample sizes sufficient
- **Diminishing returns**: Limited benefit beyond 80 clusters

### Method Performance
- Bayesian methods provide more nuanced uncertainty quantification
- Hierarchical priors improve estimation of interaction effects
- Frequentist model selection tends to be conservative

## Contact

**Keith Goldfeld**  
NYU Grossman School of Medicine  
Email: keith.goldfeld@nyulangone.org

---

*For questions about implementation or methodology, please open an issue or contact the corresponding author.*

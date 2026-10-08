# CaPri

CaPri is an R package for Bayesian Mendelian randomization with potentially pleiotropic genetic instruments.

Mendelian randomization uses genetic variants as instrumental variables to estimate the causal effect of an exposure on an outcome. However, genetic instruments may violate the exclusion restriction if they affect the outcome through pathways other than the exposure. CaPri is designed to model this situation by jointly estimating the causal effect and the direct pleiotropic effects of genetic instruments.

CaPri is fitted using GWAS summary statistics.

## Main features

CaPri estimates:

* `q_delta`: posterior probability of a causal effect;
* `hat_gamma`: estimated causal effect;
* `pip_c`: posterior probability that each instrument is pleiotropic;
* `mu_betaY`: estimated magnitude of the pleiotropic effect.

The package provides the main function:

```r
fit_capri_sumstats()
```

which fits CaPri using summary statistics.

## Installation

You can install the development version of CaPri from GitHub with:

```r
install.packages("remotes")
remotes::install_github("OuRen0104/CaPri")
```

Then load the package:

```r
library(CaPri)
```

## Real data application: Serum urate and gout

This example demonstrates how to apply CaPri to real GWAS summary statistics. We investigate the causal relationship between serum urate and gout using genetic instruments.

CaPri estimates the causal effect and pleiotropic effect of genetic instruments.

### Load real GWAS summary statistics

The `Urate_Gout` dataset contains harmonized GWAS summary statistics for serum urate (exposure) and gout (outcome).

```r
library(CaPri)
# Load data
data("Urate_Gout")

# Keep SNPs passing harmonization filters
dat <- subset(
  Urate_Gout,
  mr_keep
)

# Inspect the summary statistics
head(
  dat[, c(
    "SNP",
    "effect_allele.exposure",
    "other_allele.exposure",
    "effect_allele.outcome",
    "other_allele.outcome",
    "beta.exposure",
    "beta.outcome",
    "eaf.outcome"
  )]
)
```

The main variables used by CaPri are:

* `beta.exposure`: SNP-exposure association estimates for serum urate.
* `beta.outcome`: SNP-outcome association estimates for gout.
* `eaf.outcome`: effect allele frequency used for allele orientation.
* `n`: effective sample size of the outcome GWAS (see below).

### Prepare inputs for CaPri

CaPri requires SNP-exposure association estimates (`betaX`), SNP-outcome association estimates (`muY`), and the outcome GWAS sample size (`n`). Since gout is a binary outcome, we use the effective sample size `4 * n_case * n_control / (n_case + n_control)`, with 2,115 cases and 67,259 controls.

```r
# SNP-exposure association estimates
betaX <- as.numeric(
  dat$beta.exposure
)

# SNP-outcome association estimates
muY <- as.numeric(
  dat$beta.outcome
)

# Effect allele frequency
eaf <- as.numeric(
  dat$eaf.outcome
)

# Effective sample size of the outcome GWAS (binary outcome)
n_case <- 2115
n_control <- 67259
n <- 4 * n_case * n_control / (n_case + n_control)

n
```

### Fit CaPri

CaPri allows users to specify whether summary statistics should be used as supplied or converted to the minor allele direction.

By default, CaPri uses the supplied effect estimates without modification:
```r
allele_mode = "as_is"
```

Alternatively, users can set:
```r
allele_mode = "minor"
```

When using `allele_mode = "minor"`, effect allele frequencies (`eaf`) are required. CaPri automatically flips the signs of SNP effects when the supplied effect allele corresponds to the major allele.
For this real-data example, we use the minor allele orientation:

```r
fit <- fit_capri_sumstats(
  muY = muY,
  betaX = betaX,
  n = n,
  allele_mode = "minor",
  eaf = eaf,
  sigma_Y = 1.0,
  verbose = FALSE
)
```

### Interpret CaPri results

A fitted CaPri object contains posterior estimates of the causal effect and horizontal pleiotropic effects.

#### Posterior probability of a causal effect

```r
fit$q_delta
```

`q_delta` represents the posterior probability that a causal effect exists between serum urate and gout.
A value close to 1 indicates strong posterior evidence supporting a causal relationship.

#### Estimated causal effect

```r
fit$hat_gamma
```

`hat_gamma` is the estimated causal effect of serum urate on gout after accounting for potential horizontal pleiotropic effects.

#### Posterior probability of pleiotropic instruments

```r
fit$pip_c
```

`pip_c` gives the posterior probability that each genetic instrument has a direct pleiotropic effect on the outcome.
Genetic instruments with high `pip_c` values are more likely to violate the exclusion restriction assumption and may affect gout through pathways other than serum urate.


#### Estimated pleiotropic effect

```r
fit$mu_betaY
```

`mu_betaY` is the estimated pleiotropic effect of genetic instruments on gout.

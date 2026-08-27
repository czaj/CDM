# CDM
Count Data Models

This is a package of Matlab scripts and functions that allow for the estimation of count data models, in particular for single-site travel cost method (TCM) studies of recreation demand. The dependent variable is a count — the number of trips a respondent takes — and the package covers the sampling schemes that on-site and off-site recreation surveys typically produce.

This package includes the following models:
- Poisson regression
- Negative Binomial regression (NB2)
- Zero-inflated Poisson (ZIP)
- Zero-inflated Negative Binomial (ZINB)

Each of them can be combined with:
- zero truncation (on-site samples, where non-visitors are unobserved)
- zero truncation with endogenous stratification (on-site samples, where the probability of being sampled is proportional to the number of trips)
- right censoring at a user-specified count
- covariates of overdispersion (the NB dispersion parameter can be made a function of explanatory variables)
- weighted maximum likelihood
- analytical gradients, several Hessian and standard error options (including robust standard errors)
- equality restrictions on parameters

The models are estimated using the maximum likelihood method. Consumer surplus per trip is reported alongside the coefficients.

## Getting started

Add the `CDM` and [`Tools`](https://github.com/czaj/Tools) folders to your Matlab path. `Tools` is required — it provides the output formatting, numerical differentiation and helper functions shared with the [`DCE`](https://github.com/czaj/DCE) package.

```matlab
addpath('...\CDM')
addpath(genpath('...\Tools'))
```

`CDM_demo.m` walks through the whole family of models on simulated data and is the fastest way to see what the package does.

A minimal session looks like this:

```matlab
INPUT.Y  = trips;              % counts, one row per observation
INPUT.Xa = [ones(n,1), tc, income];   % explanatory variables, including the constant

EstimOpt.NP  = n;              % number of respondents
EstimOpt.NCT = 1;              % observations per respondent
EstimOpt.NB  = 1;              % 0 = Poisson, 1 = Negative Binomial
EstimOpt.NamesA = {'Cons';'Travel cost';'Income'};
EstimOpt.TC_par_index = 2;     % which parameter of Xa is the travel cost

[INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.NB = CDM(INPUT,Results,EstimOpt,OptimOpt);
```

`DataCleanCDM` fills in the estimation defaults, handles weights and missing data, and estimates the constants-only model used for McFadden's pseudo-R². `CDM` estimates the model you specified and prints the results (also saving them to an `.xls` file).

## Data

| field | description |
| --- | --- |
| `INPUT.Y` | counts (non-negative integers), `NCT*NP` rows |
| `INPUT.Xa` | explanatory variables of the count equation — include your own constant |
| `INPUT.Xnb` | optional covariates of overdispersion; a constant is prepended automatically |
| `INPUT.Xzinf` | covariates of the zero-inflation equation — include your own constant |
| `INPUT.W` | optional weights, one per respondent or one per observation |
| `INPUT.MissingInd` | optional 0/1 indicator; flagged observations are removed |

## Model options

Set these in `EstimOpt` before calling `DataCleanCDM`.

| option | description |
| --- | --- |
| `NP` | number of respondents |
| `NCT` | observations per respondent (1 for a cross-section) |
| `NB` | 0 = Poisson (default), 1 = Negative Binomial |
| `Zinf` | 1 for a zero-inflated model (requires `INPUT.Xzinf`) |
| `probit` | 1 to use a probit rather than a logit inflation equation |
| `Truncated` | 0 = none (default), 1 = zero truncation, 2 = zero truncation with endogenous stratification |
| `Censored` | right-censoring point; observations above it are censored (with a warning) |
| `NamesA`, `NamesNB`, `NamesZinf` | variable labels for the output |
| `TC_par_index` | index of the travel cost parameter in `Xa`, used for consumer surplus |
| `Display` | 1 (default) prints the results, 0 runs quietly |
| `ProjectName` | appended to the name of the saved `.xls` file |

## Estimation options

These are shared with the [`DCE`](https://github.com/czaj/DCE) package and behave the same way.

| option | description |
| --- | --- |
| `NumGrad` | 0 = analytical gradient (default), 1 = numerical |
| `HessEstFix` | 0 = from the optimizer (default), 1 = BHHH, 2 = high-precision BHHH, 3 = numerical, 4 = analytical (not derived for count models — falls back to 3) |
| `ApproxHess` | 1 for a BHHH approximation of the user-supplied Hessian |
| `RobustStd` | 1 for sandwich (robust) standard errors |
| `RealMin` | 1 to guard the likelihood against underflow to zero |
| `ConstVarActive`, `BActive` | impose equality restrictions on parameters |

## Results

`CDM` returns a structure with, among others:

| field | description |
| --- | --- |
| `bhat`, `std`, `R` | estimates, standard errors and p-values |
| `DetailsA` | the count equation |
| `DetailsNB` | the overdispersion equation (NB only) |
| `DetailsZinf` | the zero-inflation equation (ZIP/ZINB only) |
| `DetailsCS` | consumer surplus per trip, −1/β<sub>tc</sub>, with a delta-method standard error |
| `stats` | LL, LL of the constants-only model, McFadden's pseudo-R², AIC/n, BIC/n, n, r, k |
| `LL`, `hess`, `ihess` | log-likelihood and the (inverse) Hessian |
| `R_out` | the formatted output table, as written to `.xls` |

Note that the intercept of the overdispersion equation is reported without significance stars: H<sub>0</sub>: ln α = 0 is not the Poisson boundary, so the usual test does not answer the question people expect it to. Covariates of overdispersion are reported with ordinary tests.

## Tests

`tests/run_CDM_tests.m` runs the package test suite: analytical against numerical gradients over the whole grid of model specifications, the structural paths of `CDM` and `DataCleanCDM`, and parameter recovery on simulated data.

```matlab
run_CDM_tests
```

#
The codes are published under a Creative Commons Attribution 4.0 License. This means that you are free to use, share, or modify these codes for any purpose, even commercially. What we ask in return is that you acknowledge the source of the codes or reference one of our papers (see czaj.org/research for details).

We provide these codes for two reasons:
- Evolution - feel free to study, apply, extend, and build upon what we have accomplished.
- Efficiency - we have put considerable effort into making the codes fast and efficient. We hope to receive feedback; if you have any suggestions for improving these codes or simply making them more elegant – please let us know

The codes come with no warranty – we try to make them error free and as researcher friendly as possible; however, certain errors may remain. The demos and documentation are rather scant; if you want to use these codes, be prepared to spend a significant amount of time understanding them.

You may want to begin by adding the CDM and [Tools](https://github.com/czaj/Tools) folders to your Matlab path and checking out `CDM_demo.m`, or refer to http://czaj.org/research/supplementary-materials for software codes that accompany our papers.

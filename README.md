# CDM
Count Data Models

This is a package of Matlab scripts and functions that allow for the estimation of count data models, in particular for single-site travel cost method (TCM) studies of recreation demand. The dependent variable is a count — the number of trips a respondent takes — and the package covers the sampling schemes that on-site and off-site recreation surveys typically produce.

This package includes the following models:
- Poisson regression
- Negative Binomial regression: NB2 (variance λ + αλ²), NB1 (variance λ(1+α)) and NB-P (variance λ + αλ^P, with the power estimated)
- Zero-inflated Poisson (ZIP) and Negative Binomial (ZINB)
- Hurdle models: a binary participation equation plus a zero-truncated count

Each of them can be combined with:
- zero truncation (on-site samples, where non-visitors are unobserved)
- zero truncation with endogenous stratification (on-site samples, where the probability of being sampled is proportional to the number of trips)
- right censoring at a user-specified count
- covariates of overdispersion (the NB dispersion parameter can be made a function of explanatory variables)
- weighted maximum likelihood
- analytical gradients, several Hessian and standard error options (including robust standard errors)
- equality restrictions on parameters

The models are estimated using the maximum likelihood method, with analytical gradients throughout. Consumer surplus per trip is reported alongside the coefficients.

All three Negative Binomial parameterisations are the same likelihood with θ = λ<sup>2−P</sup>/α: NB2 is P = 2, NB1 is P = 1, and NB-P estimates P. Since NB-P nests the other two, it is the honest way to choose between them — but the power is not always well identified, so check its standard error.

## Getting started

Add the `CDM` and [`Tools`](https://github.com/czaj/Tools) folders to your Matlab path. **`Tools` is a hard dependency** — CDM carries no output code of its own, so the printed tables, the `.xls` export, the numerical differentiation and the Hessian diagnostics all come from `Tools`, exactly as in the [`DCE`](https://github.com/czaj/DCE) package.

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

`DataCleanCDM` sets up an `output/` folder next to the calling script, starts a diary of the session there, fills in the estimation defaults, handles weights and missing data, and estimates the constants-only model (returned as `Results.CDM0`) whose log-likelihood gives McFadden's pseudo-R². `CDM` estimates the model you specified, prints the results and saves them to `output/` as an `.xlsx` file named after the model — `NB_results_demo.xlsx`, `ZINB_results_demo.xlsx` and so on. All of that comes from `Tools`, shared with the `DCE` package.

Note that `DataCleanCDM` returns a fresh `Results`, so when you estimate several models in one script keep your own accumulator:

```matlab
[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0;
Results.NB = CDM(INPUT,Results,EstimOpt,OptimOpt);
```

## Data

| field | description |
| --- | --- |
| `INPUT.Y` | counts (non-negative integers), `NCT*NP` rows |
| `INPUT.Xa` | explanatory variables of the count equation — include your own constant |
| `INPUT.Xnb` | optional covariates of overdispersion; a constant is prepended automatically |
| `INPUT.Xzinf` | optional covariates of the zero-inflation equation; a constant is prepended automatically |
| `INPUT.W` | optional weights, one per respondent or one per observation |
| `INPUT.MissingInd` | optional 0/1 indicator; flagged observations are removed |

## Model options

Set these in `EstimOpt` before calling `DataCleanCDM`.

| option | description |
| --- | --- |
| `NP` | number of respondents |
| `NCT` | observations per respondent (1 for a cross-section) |
| `NB` | 0 = Poisson (default), 1 = Negative Binomial |
| `NBtype` | 2 = NB2, variance λ + αλ² (default); 1 = NB1, variance λ(1+α); 3 = NB-P, variance λ + αλ^P with `P` estimated |
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
| `stats` | LL, LL of the constants-only model, McFadden's pseudo-R², (unused), AIC/n, BIC/n, n, r, k |
| `LLdetailed` | weighted log-likelihood per observation |
| `HessDiagnostics` | rank, condition number and inversion method of the Hessian |
| `LL`, `hess`, `ihess` | log-likelihood and the (inverse) Hessian |
| `R_out` | the formatted output table, as written to `.xls` |

Note that the intercept of the overdispersion equation is reported without significance stars: H<sub>0</sub>: ln α = 0 is not the Poisson boundary, so the usual test does not answer the question people expect it to. Covariates of overdispersion are reported with ordinary tests.

## Hurdle models

`CDM_hurdle` estimates a two-part model: a binary equation for whether the respondent participates at all, and a zero-truncated count for how many trips a participant takes.

```matlab
INPUT.Xh = [];               % participation covariates; defaults to Xa
EstimOpt.NamesH = {'Travel cost';'Income'};
Results.Hurdle = CDM_hurdle(INPUT,Results,EstimOpt,OptimOpt);
```

For recreation demand this is usually better motivated than zero inflation: not going at all is a different decision from deciding how often to go, rather than a hidden class. The log-likelihood separates,

```
LL = sum_{y=0} log(1-P) + sum_{y>0} [ log(P) + log f(y | y>0) ]
```

and the two parts share no parameters, so each is maximised on its own and the information matrix is block diagonal. The result is shaped like any other CDM result — the participation equation is stored in the zero-inflation slot, so `CDM_post` and `CDM_compare` work on it unchanged, and expected trips come out as P·E[Y|Y>0].

## Travel cost

`TravelCost` builds the travel cost variable. The money part is uncontroversial; the value of travel time is not, which is why the rule is an explicit choice and why the default returns a whole sensitivity panel rather than one number.

```matlab
Trip = struct('Distance',dist,'CostPerKm',0.35,'Time',dist/70,'Wage',wage);
[TC,D] = TravelCost(Trip);              % columns for 0, 1/4, 1/3 and 1/2 of the wage
INPUT.Xa = [ones(n,1), TC(:,3), inc];   % 1/3 of the wage
```

| `Opt.TimeRule` | what it does |
| --- | --- |
| `'fraction'` | a fraction of the wage — one column per value in `Opt.TimeFraction`, default `[0 1/4 1/3 1/2]` |
| `'full'` | the full wage rate |
| `'full_or_zero'` | the full wage for respondents who would have preferred a shorter trip, nothing for those who would not. Needs `Trip.ValuesTime`. Travel that someone would not shorten is not time given up, so charging them for it overstates the price they faced |
| `'none'` | money cost only |

Distance and time are doubled unless `Trip.RoundTrip = 1`. `Trip.OnSiteTime` is valued by the same rule.

## Post-estimation

`CDM_post` turns a fitted model into the numbers a travel cost study actually reports, and prints them:

```matlab
Results.NB = CDM_post(Results.NB);
```

| what | detail |
| --- | --- |
| `EY` | expected trips per observation, matched to the sampling scheme — `lambda` for a general population sample, `lambda/(1-P(0))` under zero truncation, `E[Y^2]/E[Y]` under endogenous stratification (see `CDM_EY`) |
| `CS.PerTrip` | consumer surplus per trip, −1/β<sub>tc</sub>, with a delta-method standard error, a Krinsky-Robb interval and the interval obtained by inverting the confidence interval on β (the exact interval for a ratio with a constant numerator) |
| `CS.PerPerson` | seasonal surplus per person, E[Y]·(−1/β<sub>tc</sub>), with a Krinsky-Robb interval |
| `Margins` | average marginal effects on E[Y] and elasticities, with Krinsky-Robb standard errors. Elasticities are reported only for strictly positive regressors |
| `CountFit` | observed against predicted count frequencies and a Pearson chi-square. Predicted probabilities come from the estimated likelihood, so truncation, stratification, censoring and zero inflation are all handled consistently |

Right censoring deliberately does not enter `E[Y]`: it describes how trips were recorded, not how many were taken, and welfare measures need the latent expectation.

`CDM_compare` tests one model against another and picks the right test:

```matlab
CDM_compare(Results.POISS,Results.NB);   % likelihood ratio, boundary-corrected
CDM_compare(Results.POISS,Results.ZIP);  % Voung, non-nested
```

Poisson is Negative Binomial with α = 0, which sits on the boundary of the parameter space, so the likelihood ratio statistic follows an equal mixture of a point mass at zero and χ²(1) — using plain χ²(1) makes the test about twice as conservative as it should be. Zero inflation is not identified under its own null, so `CDM_compare` falls back to the Voung test there.

## Tests

`tests/run_CDM_tests.m` runs the package test suite: analytical against numerical gradients over the whole grid of model specifications, the structural paths of `CDM` and `DataCleanCDM`, parameter recovery on simulated data, and the post-estimation results — including a check of every expected-trip formula against Monte Carlo.

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

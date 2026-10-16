# Hybrid travel-cost/count models

`CDM_hybrid` jointly estimates a count equation and an ordinal measurement
model with normally distributed latent factors. This extension does not change
the existing `CDM` and `DataCleanCDM` estimators.

## Model

For respondent i and factor l:

```
LV_il = Xstr_i * gamma_l + epsilon_il,     epsilon_i ~ N(0, I)
log(mu_i) = Xa_i * beta + sum_l LV_il * (Xlatent_i * delta_l)
Y_i | LV_i ~ Poisson(mu_i) or NB2(mu_i, alpha)
P(I_ij <= k | LV_i) = logistic(tau_jk - loading_j * LV_i,factor(j))
```

The model integrates the product of the count likelihood and all observed
indicator likelihoods. It never treats estimated factor scores as error-free
count regressors. Fixed quadrature makes the objective and score deterministic.
Gauss-Hermite is the default; a uniform midpoint grid with normal-density
weights is available for sharply concentrated indicator likelihoods.
Missing indicators contribute one to the likelihood, under an
ignorable missingness assumption; missing count/design values must be addressed
before estimation.

Identification uses unit residual variance and no structural intercept. Each
item loads on one factor, each factor has at least three indicators, and its
first loading is constrained positive to fix orientation. Ordered thresholds
replace item intercepts. Independent latent disturbances are a maintained
restriction, not evidence that the underlying motivations are independent.
Three indicators alone do not establish construct validity or empirical
identification: inspect loadings, fit, Hessian diagnostics and quadrature
sensitivity.

## Sampling matters

With `Truncated=1` and the default `LatentDistribution='population'`, the
sample is selected on **this outcome** being positive:

```
L_i = integral f(Y_i | LV) p(I_i | LV) phi(LV | Xstr_i) dLV
      / integral [1 - f(0 | LV)] phi(LV | Xstr_i) dLV.
```

The denominator is outside the integral. Integrating a conditionally truncated
count likelihood against the population factor distribution would ignore how
selection changes the distribution of motivations. The implemented correction
does not correct selection on a different recall period, panel attrition, quota
sampling or on-site length bias. `Truncated=2`, panel records, censoring,
zero inflation and NB1/NB-P are not supported by this extension.

An explicit sensitivity, `LatentDistribution='selected'`, instead assumes
the normal factor law **among visitors**, including its structural equation:

```
L_i = integral [f(Y_i | LV) / (1-f(0 | LV))] p(I_i | LV)
               phi_selected(LV | Xstr_i) dLV.
```

This is a different sampling/distribution assumption, not a numerical
approximation to the default likelihood. It requires `Truncated=1`.
Here the indicator marginal integrates just its measurement model under
the visitor-normal law; the count component does not tilt that marginal.
In the log-series boundary the conditional kernel is `p^Y/(Y*softplus(v))`.
Results carry the suffix `_SELECTED_LATENT`. Selected means and item-only
posterior scores use the visitor law; population means and participation
probabilities are not reported by this implementation. Do not apply ordinary
population-welfare formulas to these outputs. The two latent-distribution
assumptions coincide at a zero count-factor link, but their full models
are not nested and cannot be compared with a regular chi-squared LR test.

## Package layout and the standard calling script

Hybrid files follow CDM's flat layout beside the existing estimators:

- `DataCleanCDM_hybrid.m`: validates inputs, resolves options and starting values;
- `CDM_hybrid.m`: estimation, diagnostics and optional output;
- `LL_CDM_hybrid.m`: respondent log likelihood and analytical scores;
- `LL_CDM_hybrid_MATlike.m`: summed negative objective and analytical gradient;
- `CDM_hybrid_output.m`: native CDM/DCE tables and XLSX export;
- `CDM_demo_hybrid.m`: complete calling-script example.

Standard calls use the same four objects as CDM:

```matlab
addpath('.../CDM','.../Tools')
[INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM_hybrid(INPUT,EstimOpt);
Results = CDM_hybrid(INPUT,Results,EstimOpt,OptimOpt);
```

`SaveXLS`, `SaveTxtOutput`, `OutputDir`, `ProjectName`, `OutputFilePrefix`,
`SourceScript` and `OutputLogFile` control output. Explicit standalone export
is available as `Results = CDM_hybrid_output(Results,exportOpt)`, without a
refit. `Results.OutputXLSX` records the actual file. Export uses shared
`Tools/genOutput` and `writeResultsXlsx`; the `Results`, `Word 2dp`, and
`Word 4dp` sheets retain the existing parameter/significance/SE/p format.
`Script` and `Log` are included when the corresponding files exist.

Natural dispersion, positive anchor loadings and ordered thresholds have
delta-method SEs using the full covariance matrix. A separate raw block
retains every optimizer parameter. Alpha and the constrained positive anchor
have no ordinary zero-test p-value; raw parameters have no significance stars.
Public inference requires valid convergence/Hessian diagnostics and a passed
quadrature check; unchecked local SEs are labelled diagnostic. The Hessian is
a numerical derivative of the analytical score, not an analytical Hessian.

## Inputs and a minimal call

```matlab
addpath('.../CDM')
INPUT.Y = trips;                         % n x 1 nonnegative integers
INPUT.Xa = [ones(n,1), travelCost, age];  % n x p, own intercept
INPUT.Xstr = age;                        % n x q, NO intercept
INPUT.Xmea = indicators;                 % n x J, 1..K; missing = NaN
INPUT.Xlatent = ones(n,1);               % optional, latent count intercept
EstimOpt = struct('NB',1,'Truncated',1,'MeaMatrix',ones(3,1));
Results = CDM_hybrid(INPUT,EstimOpt);
```

Raw 0..10 importance ratings must be recoded to 1..11 **after** removing any
special missing codes. For example, a raw 11 meaning "don't know" is missing,
not the highest importance category. `MeaMatrix` is J x L with one 1 per row.
`Xlatent` can contain a constant and observed moderators. Include their main
effects in `Xa` when fitting interactions; the package cannot establish a
causal interpretation from a supplied matrix.

Every fitted item must have at least two observed, consecutive categories
starting at 1. Empty categories cause collapsed/unidentified threshold gaps;
merge/recode represented ordinal levels deliberately and retain their mapping.

| Option | Default / meaning |
| --- | --- |
| `NNodes` | 15 points per factor; Hermite range 3..512; normal-grid up to the total node limit; check a higher value |
| `Quadrature` | `'hermite'` or `'normal-grid'`; fixed nodes throughout optimization |
| `GridLimit` | 12 residual standard deviations in either direction for normal-grid; check both range and resolution |
| `MaxQuadratureNodes` | 4096 total tensor points, at most three factors |
| `LogSeries` | false; selected NB2 boundary, requires `NB=0, Truncated=1`; see below |
| `LatentDistribution` | `'population'` by default; `'selected'` describes normally distributed factors among visitors and requires `Truncated=1` |
| `LatentActive` | true; false drops only the latent effects on counts |
| `Start` | optional raw starting vector for the resolved parameter layout |
| `EvaluateOnly` | false; true resolves the layout and evaluates without fitting |
| `ComputeHessian` | true; central differences of the analytical score |
| `RobustStd` | false; true uses the observation-score sandwich |
| `NamesA`, `NamesStr`, `NamesLatent`, `NamesMea` | labels of the corresponding matrices |

`Results.bhat`, `std`, `ParameterNames` and `Layout` describe raw parameters;
`Parameters` contains transformed loadings, thresholds and dispersion.
`HessDiagnostics.ValidInference` gates standard errors using convergence,
gradient, positive-definiteness and conditioning checks. `LLdetailed` and
`Score` are weighted respondent contributions. The saved `INPUT` and `EstimOpt`
can be passed to `LL_CDM_hybrid` to evaluate the log likelihood and score.
This function returns log likelihood; the older `LL_CDM` returns its negative.

Optimizer convergence and a positive Hessian do **not** establish numerical
integration accuracy. `QuadratureChecked=false` makes this distinction explicit.
Large item loadings can require far more than 15 nodes. Evaluate the same fitted
parameters at increased node counts, then refit both restricted and unrestricted
models at that accuracy and compare log likelihoods and coefficients. Logarithmic
Hermite weights preserve tails at high quadrature orders. The BlueAdapt example
implements these checks and withholds LR p-values and welfare when unverified.
For example, use `Quadrature='normal-grid', NNodes=128`, then refit at 256
with the same `GridLimit`. The grid weights are deliberately not renormalized:
`EstimOpt.QuadratureMass` exposes retained normal mass. A coarse grid or short
range can be inaccurate. Check mean **absolute** per-person LL differences
to avoid cancellation, and compare all coefficients as well as the latent
effect. These checks concern numerical accuracy, not validity of the
measurement construct or the sampling assumptions.

`EYPopulation` integrates population motivations; `EYSelected` conditions on
Y>0. `Prediction.EYGivenIndicatorsPopulation` and
`Prediction.EYGivenIndicatorsSelected` additionally condition on indicators.
`PosteriorMean` uses the observed count too. These are different quantities.
`Prediction.LogIndicatorMarginalSelected` is the marginal indicator log
likelihood among positive-count visitors. For `Truncated=1` and unit weights,
subtract it from `LLdetailed` to score the count conditional on the indicators. With weights,
first remove the observation weight from `LLdetailed`. This decomposition
helps separate improvement in the indicator distribution from improvement
in predicting trips.
The standard second argument supplies previous results; starts are matched
by parameter name. The older `CDM_hybrid(INPUT,EstimOpt,OptimOpt,previousFit)`
order and two/three-argument calls remain supported. Explicit `Start` vectors
must match the resolved layout. Passing an earlier result does not reuse its
likelihood on new data: any likelihood comparison evaluates the current input.

Set `LatentActive=false` to estimate the same joint model with the count-factor
link constrained to zero. This is the appropriate restricted model for testing
the added count-factor link: it contains the **same indicators and people**.
Do not compare a count-only likelihood/AIC with a joint count-plus-indicator
likelihood/AIC. Weighted pseudo-likelihoods require appropriate survey
inference; an ordinary likelihood-ratio p-value is not automatically valid.

Use `CDM_demo_hybrid` for a reproducible example and `tests/run_CDM_tests` for
the numerical tests. The engine requires MATLAB and Optimization Toolbox.
Predictions and latent posterior scores are diagnostic model outputs, not new
observations. Scores using Y must not be reused to claim out-of-sample prediction
or pretreatment moderation of Y.

## The selected NB2 boundary

A positive-count NB2 fit can approach `alpha -> infinity` and population
`mu -> 0`, keeping `alpha*mu` finite. Its intercept and log dispersion then
move in opposite directions along a nearly flat likelihood. A positive
optimizer exit flag does not make their separate estimates informative.
`LogSeries=true` estimates this limiting selected family directly, without
an arbitrary dispersion cap. Write `v = Xa*beta + LV*(Xlatent*delta)`,
`p = logistic(v)` and `s = log(1+exp(v))`. The joint likelihood is

```
L_i = integral [p(LV)^Y_i / Y_i] p(I_i | LV) phi(LV | Xstr_i) dLV
      / integral s(LV) phi(LV | Xstr_i) dLV,       Y_i >= 1.
```

Thus `v` represents `log(alpha*mu)` in the limit. There is no dispersion
parameter. The normal latent density is a reference distribution; the
selected latent distribution is tilted by `s(LV)`. This is not a mixture
of individually normalized log-series distributions under an unchanged
normal latent density. For a warm start from NB2, map shared parameter names
and add its log dispersion to the count intercept (when that column is one).

`EYSelected = integral exp(v) phi / integral s phi`; indicator-conditioned
predictions include measurement probabilities in both integrals. Population
means and participation probabilities are unidentified and return `NaN`.
`Prediction.SelectionWeightIntegral` is an unnormalized weight integral,
**not** a probability. Population welfare cannot be recovered from this
boundary model with the usual `-1/beta_cost` calculation.

The restricted-versus-full latent-link test compares two log-series fits
with identical indicators and respondents. An ordinary chi-squared test
for NB2 versus its boundary is not valid; report that comparison descriptively.
The tests check normalization, analytic gradients, convergence of NB2 to
this limit, and selected predictions against independent calculations.

Passing Hessian and quadrature checks establishes only local numerical
inference under the chosen model. It does not establish predictive
calibration. In particular, a heavy-tailed selected model can improve log
likelihood while substantially overpredicting the mean. Compare predicted
and observed counts, inspect the largest predictions, and evaluate held-out
respondents using predictions conditioned on indicators, never on their Y.

## Interpretation for recreation demand

This is a single-site demand model: Y counts trips to a defined destination.
It can explain frequency using latent motivations. To explain *which* beach
was chosen, use an alternative-specific location model with the available
sites, travel costs and attributes, or develop a joint frequency/site-choice
model. Do not relabel the count model as spatial choice.

Treatment-by-latent interactions belong in `Xlatent` (and treatment main effects
in `Xa`). Causal moderation requires motivations measured before randomization
and an outcome observed after treatment. Post-treatment knowledge can be a
mediator and needs a different estimand/design. Do not treat a baseline
retrospective trip count as the effect of information supplied during the
survey. If travel cost changes latent effects or the structural distribution,
the ordinary -1/beta cost welfare formula no longer suffices; this extension
deliberately does not call `CDM_post` or report automatic consumer surplus.

## Methodological background

The integrated structural/measurement/behavioral likelihood follows the
[Ben-Akiva et al. integrated choice and latent-variable framework](https://www.joanwalker.com/uploads/3/6/9/5/3695513/benakivawalkeretal_iclv_chapter_2002.pdf),
with a count outcome replacing the discrete choice component. Latent variables
have also been used in travel-cost research, including
[Englin and Shonkwiler's travel-price model (1995)](https://doi.org/10.1006/jeem.1995.1053).
That paper concerns latent travel costs; it is context, not a claim that this
implementation or the BlueAdapt motivation specification replicates it.

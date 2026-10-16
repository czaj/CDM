function Results = CDM_hybrid(INPUT,varargin)
% CDM_hybrid Joint latent-variable Poisson/NB2 travel-demand model.
%   R = CDM_hybrid(INPUT,EstimOpt,OptimOpt) estimates counts and ordered-logit
%   indicators jointly; OptimOpt is optional. The standard CDM calling form is
%   R = CDM_hybrid(INPUT,Results_old,EstimOpt,OptimOpt), using the four outputs
%   of DataCleanCDM_hybrid. Raw matrices are also accepted directly.
%   The older (INPUT,EstimOpt,OptimOpt,Results_old) order remains supported.
%   Results_old may be a previous hybrid fit or a structure containing fits.
%
%   INPUT.Y       n-by-1 nonnegative integer trip counts.
%   INPUT.Xa      n-by-k count covariates, INCLUDING any desired intercept.
%   INPUT.Xmea    n-by-J ordinal indicators, coded consecutively from 1;
%                 NaN indicators are integrated out (ignorable missingness).
%   INPUT.Xstr    n-by-q latent-mean covariates, WITHOUT any constant; optional.
%   INPUT.Xlatent n-by-m modifiers of latent count effects; default ones(n,1).
%                 Include a constant here to estimate a main latent effect.
%   INPUT.W       optional positive observation weights.
%
%   LV = Xstr*Gamma + epsilon, epsilon ~ N(0,I).
%   log E[Y|LV] = Xa*beta + sum_l LV_l .* (Xlatent*delta_l).
%   P(I_j<=k|LV) = logistic(tau_jk - loading_j*LV_factor(j)).
%   Indicators are conditionally independent given LV. The first loading of
%   each factor is exp(raw loading), fixing orientation; other loadings are
%   unconstrained. The first threshold is free and later gaps are exp(raw gap).
%
%   EstimOpt.MeaMatrix: J-by-L binary simple structure, one factor per item,
%       at least three informative indicators per factor (default one factor).
%   NB=0/1 (Poisson/NB2); NBtype=2 only. Truncated=0/1. For Truncated=1,
%       sampling is random from people with Y>0. The JOINT numerator is
%       divided by the population-integrated P(Y>0|covariates). This is not
%       a mixture of conditionally truncated count distributions.
%   LatentActive=false removes delta but retains the same measurement model,
%       enabling a nested joint-likelihood test of the latent count effects.
%   NNodes=15 nodes per dimension (3..512), at most three independent factors.
%       Increase and refit until likelihood and parameters are stable; a
%       positive Hessian alone does not establish quadrature accuracy.
%   Quadrature='hermite' (default) or 'normal-grid'. The latter uses midpoint
%       nodes on [-GridLimit,GridLimit], GridLimit=12, with normal-density
%       integration weights and no finite-range renormalization. Its NNodes
%       limit is MaxQuadratureNodes (4096 by default); check resolution AND
%       the integration range. Nodes remain fixed during optimization.
%   LogSeries=true fits the selected NB2 boundary as alpha->Inf at fixed
%       alpha*mu; requires NB=0 and Truncated=1. The count index then means
%       log(alpha*mu), not log(population mean). The joint numerator kernel is
%       p^Y/Y and the denominator integrates -log(1-p), p=logistic(index).
%       Population mean and participation probabilities are not identified
%       by this limiting selected model and are reported as NaN.
%   LatentDistribution='population' (default) places the normal factor law
%       before positive-count sampling and normalizes outside the integral.
%       'selected' instead assumes normal factors among visitors and mixes
%       the positive-count conditional PMF inside the integral. This is a
%       different sampling assumption, not a correction to the default.
%       It requires Truncated=1; population means/probabilities remain NaN.
%   NCategories optionally fixes the known number of categories per item;
%       all categories must occur to avoid unidentified/collapsed thresholds.
%   Start optionally supplies the full raw parameter vector. EvaluateOnly=true
%       resolves inputs/layout and evaluates Start without optimization.
%   SaveXLS and SaveTxtOutput default to Display~=0, using the companion Tools
%       output routines. ProjectName, OutputDir, SourceScript, OutputLogFile
%       and OutputFilePrefix control the standard output. Set both save flags
%       to 0 for an estimation-only run. EvaluateOnly never exports files.
%   ComputeHessian=true computes an observed numerical Hessian of the analytic
%       score. RobustStd=false; if true use observation-score sandwich SEs.
%       Weights enter the objective as supplied; sandwich uses squared weights.
%
%   R.bhat and R.std refer to RAW parameters, labelled in R.ParameterNames;
%   R.Parameters contains transformed loadings, thresholds and dispersion.
%   R.PosteriorMean conditions on counts AND indicators. R.EYPopulation and
%   R.EYSelected integrate the specified latent distribution and do NOT use
%   posterior scores. See R.Prediction for indicator-conditioned predictions.
%   Welfare/consumer surplus is intentionally not reported: interactions and
%   cost-dependent latent means require an explicit welfare specification.
%
%   No panel, censoring, endogenous stratification, zero inflation, factor
%   correlations, random loading scales or NB1/NB-P are silently approximated.
%   See also DataCleanCDM_hybrid, LL_CDM_hybrid, LL_CDM_hybrid_MATlike,
%       CDM_hybrid_output, CDM.

[EstimOpt,OptimOpt,Results_old] = argumentsForHybrid(varargin{:});
[INPUT,~,EstimOpt,DefaultOptimOpt] = DataCleanCDM_hybrid(INPUT,EstimOpt);
estimationClock = tic;
if isfield(EstimOpt,'Start') && ~isempty(EstimOpt.Start)
    b0 = EstimOpt.Start(:);
else
    b0 = previousStart(Results_old,INPUT,EstimOpt);
end
if numel(b0) ~= EstimOpt.Params || any(~isfinite(b0)) || ~isreal(b0)
    error('CDM_hybrid:Start','Starting values must be a finite real vector of length %d.',EstimOpt.Params);
end

if EstimOpt.EvaluateOnly
    bhat = b0;
    exitflag = 0;
    output = struct('message','EvaluateOnly: no estimation or inference performed.');
    OptimOptUsed = OptimOpt;
else
    if isempty(OptimOpt)
        OptimOptUsed = DefaultOptimOpt;
    elseif isstruct(OptimOpt)
        OptimOptUsed = optimset(OptimOpt,'GradObj','on');
    else
        OptimOptUsed = optimoptions(OptimOpt,'SpecifyObjectiveGradient',true);
    end
    [f0,~] = LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOptUsed,b0);
    if ~isfinite(f0)
        error('CDM_hybrid:Start','Starting values produce a nonfinite likelihood. Rescale covariates or change Start.');
    end
    [bhat,~,exitflag,output] = fminunc(@(b) LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOptUsed,b),b0,OptimOptUsed);
end

[lli,score,detail] = LL_CDM_hybrid(bhat,INPUT,EstimOpt);
p = numel(bhat);
H = nan(p);
V = nan(p);
gradient = -sum(score,1)';
diagnostics = struct('Method','central differences of analytic score', ...
    'PositiveDefinite',false,'ReciprocalCondition',NaN, ...
    'Converged',exitflag>0,'GradientPerWeight',norm(gradient,Inf)/sum(INPUT.W), ...
    'ValidInference',false);
if EstimOpt.ComputeHessian && ~EstimOpt.EvaluateOnly && all(isfinite(lli))
    H = zeros(p);
    for j = 1:p
        step = eps^(1/3)*(1+abs(bhat(j)));
        bp = bhat; bm = bhat;
        bp(j) = bp(j)+step; bm(j) = bm(j)-step;
        [~,gp] = LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOptUsed,bp);
        [~,gm] = LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOptUsed,bm);
        H(:,j) = (gp-gm)/(2*step);
    end
    H = (H+H')/2;
    if all(isfinite(H(:)))
        [~,cholflag] = chol(H);
        diagnostics.PositiveDefinite = cholflag == 0;
        diagnostics.ReciprocalCondition = rcond(H);
        diagnostics.ValidInference = diagnostics.PositiveDefinite && ...
            diagnostics.ReciprocalCondition > 1e-10 && ...
            diagnostics.Converged && diagnostics.GradientPerWeight < 1e-3;
        if diagnostics.ValidInference
            V = H\eye(p);
            if EstimOpt.RobustStd, V = V*(score'*score)*V; end
            V = (V+V')/2;
        end
    end
    if ~diagnostics.ValidInference
        warning('CDM_hybrid:Inference','SEs withheld: optimizer convergence, score or Hessian checks failed. Inspect HessDiagnostics.');
    end
end

Results = struct();
Results.Model = 'HYBRID_POISSON';
if EstimOpt.NB, Results.Model = 'HYBRID_NB2'; end
if EstimOpt.LogSeries, Results.Model = 'HYBRID_LOGSERIES_BOUNDARY'; end
if strcmp(EstimOpt.LatentDistribution,'selected'), Results.Model = [Results.Model '_SELECTED_LATENT']; end
Results.bhat = bhat;
Results.Start = b0;
Results.std = sqrt(diag(V));
Results.ihess = V;
Results.hess = H;
Results.HessDiagnostics = diagnostics;
Results.QuadratureChecked = false;
Results.R = [bhat,Results.std,erfc(abs(bhat./Results.std)/sqrt(2))];
Results.LL = sum(lli);
Results.LLdetailed = lli;
Results.Score = score;
Results.gradient = gradient;
Results.AIC = 2*p-2*Results.LL;
Results.BIC = log(numel(INPUT.Y))*p-2*Results.LL;
Results.exitflag = exitflag;
Results.exitf = exitflag;
Results.output = output;
Results.ParameterNames = EstimOpt.ParameterNames;
Results.Layout = EstimOpt.Layout;
Results.Parameters = detail.Parameters;
Results.PosteriorMean = detail.PosteriorMean;
Results.PosteriorMeanIndicators = detail.PosteriorMeanIndicators;
Results.EYPopulation = detail.EYPopulation;
Results.EYSelected = detail.EYSelected;
Results.EY = detail.EYPopulation;
if EstimOpt.Truncated, Results.EY = detail.EYSelected; end
Results.Prediction = detail;
Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOptUsed;
Results.output_txt = '';
if isfield(EstimOpt,'OutputLogFile'), Results.output_txt = EstimOpt.OutputLogFile; end
Results.InferenceNote = ['AIC/BIC use the joint count-indicator likelihood; compare the same indicators and observations. ' ...
    'BIC uses the number of independent rows. Arbitrarily rescaled survey weights affect these criteria. ' ...
    'SEs are on the raw parameter scale and do not include quadrature approximation error. ' ...
    'Quadrature accuracy has not been checked by this fit; increase NNodes and refit before interpreting inference.'];
if EstimOpt.LogSeries
    Results.InferenceNote = [Results.InferenceNote ' The log-series family is a selected NB2 boundary; ' ...
        'an ordinary NB2-versus-boundary likelihood-ratio reference distribution is not valid. ' ...
        'The latent count-link test must compare restricted/full fits within the same family.'];
end
if strcmp(EstimOpt.LatentDistribution,'selected')
    Results.InferenceNote = [Results.InferenceNote ' The normal latent distribution and structural equation refer to visitors. ' ...
        'Population-normal and selected-normal alternatives are not nested; do not use a regular LR test to compare them.'];
end
if EstimOpt.Display
    fprintf('%s: n=%d, factors=%d, quadrature=%d, parameters=%d\n', ...
        Results.Model,numel(INPUT.Y),EstimOpt.NLatent,numel(EstimOpt.QuadratureWeights),p);
    fprintf('Joint LL %.6f; AIC %.3f; BIC %.3f; exitflag %d\n',Results.LL,Results.AIC,Results.BIC,exitflag);
end
Results.EstimationSeconds = toc(estimationClock);
if EstimOpt.EvaluateOnly, Results.EstimationSeconds = NaN; end
if ~EstimOpt.EvaluateOnly && EstimOpt.SaveXLS
    Results = CDM_hybrid_output(Results);
end
end

function [O,Optim,R] = argumentsForHybrid(varargin)
O = struct(); Optim = []; R = [];
if numel(varargin)>3
    error('CDM_hybrid:Arguments','Use INPUT, EstimOpt, optional OptimOpt, or the standard four-argument CDM form.');
end
if numel(varargin)<3
    if ~isempty(varargin), O = varargin{1}; end
    if numel(varargin)==2, Optim = varargin{2}; end
else
    a = varargin{1}; b = varargin{2}; c = varargin{3};
    % An optimizer in argument 3 identifies the legacy order; an optimizer
    % in argument 4 identifies the standard CDM order. Numeric Display is a
    % model option; textual Display ('off', 'iter', ...) is an optimizer option.
    legacy = isResult(c) || isOptimizer(b);
    if ~legacy && ~isResult(a) && ~isempty(a) && isstruct(a) && ...
            ~isempty(fieldnames(a)) && ~isOptimizer(c)
        legacy = isEstimation(a) && ~isEstimation(b);
    end
    if legacy, O = a; Optim = b; R = c;
    else, R = a; O = b; Optim = c;
    end
end
if isempty(O), O = struct(); end
end

function yes = isOptimizer(value)
yes = ~isempty(value) && ~isstruct(value);
if ~isstruct(value) || isempty(value), return; end
fields = {'Algorithm','GradObj','SpecifyObjectiveGradient','MaxIter', ...
    'MaxIterations','MaxFunEvals','MaxFunctionEvaluations','TolFun','TolX', ...
    'OptimalityTolerance','StepTolerance','Hessian','HessianApproximation'};
yes = any(isfield(value,fields)) || (isfield(value,'Display') && ...
    (ischar(value.Display) || isstring(value.Display)));
end

function yes = isEstimation(value)
yes = isstruct(value) && any(isfield(value,{'NB','Truncated','LatentActive', ...
    'NNodes','Quadrature','Display','EvaluateOnly','Start','ProjectName', ...
    'MeaMatrix','LogSeries','LatentDistribution','ComputeHessian','SaveXLS', ...
    'SaveTxtOutput','NamesA','NamesStr','NamesLatent','NamesMea','NCategories', ...
    'NBtype','NCT','NLatent','Censored','Zinf','RobustStd','GridLimit', ...
    'MaxQuadratureNodes','OutputDir','OutputFilePrefix','SourceScript'}));
end

function yes = isResult(value)
yes = false;
if ~isstruct(value) || ~isscalar(value), return; end
if isfield(value,'bhat'), yes = true; return; end
fields = fieldnames(value);
for j = 1:numel(fields)
    fit = value.(fields{j});
    if isstruct(fit) && isscalar(fit) && isfield(fit,'bhat'), yes = true; return; end
end
end

function b = previousStart(R,I,O)
b = O.InitialValues;
if isempty(R) || (isstruct(R) && isscalar(R) && isempty(fieldnames(R))), return; end
if ~isstruct(R) || ~isscalar(R)
    error('CDM_hybrid:Start','Results_old must be a previous hybrid fit or a scalar results structure.');
end
if isfield(R,'bhat')
    b = mappedStart(R,b,O);
    return;
end
% A standard Results accumulator can hold restricted and joint fits. Compare
% their mapped starts at the CURRENT likelihood, rather than stored LLs that
% may refer to a different sample, quadrature rule or parameter dimension.
bestLL = sum(LL_CDM_hybrid(b,I,O));
fields = fieldnames(R);
for j = 1:numel(fields)
    fit = R.(fields{j});
    if ~isstruct(fit) || ~isscalar(fit) || ~isfield(fit,'bhat') || ~sameFamily(fit,O), continue; end
    candidate = mappedStart(fit,b,O);
    if numel(candidate)~=O.Params || any(~isfinite(candidate)) || ~isreal(candidate), continue; end
    candidateLL = sum(LL_CDM_hybrid(candidate,I,O));
    if isfinite(candidateLL) && candidateLL>bestLL
        b = candidate; bestLL = candidateLL;
    end
end
end

function b = mappedStart(R,b,O)
if ~sameFamily(R,O)
    error('CDM_hybrid:Start','Previous fit uses a different count family, truncation or latent-distribution assumption. Supply an explicit transformed Start instead.');
end
if isfield(R,'ParameterNames')
    if numel(R.ParameterNames)~=numel(R.bhat)
        error('CDM_hybrid:Start','Previous fit has inconsistent parameter names and bhat.');
    end
    if isequal(O.ParameterNames,R.ParameterNames)
        b = R.bhat(:); return;
    end
    if numel(unique(O.ParameterNames))~=numel(O.ParameterNames) || ...
            numel(unique(R.ParameterNames))~=numel(R.ParameterNames)
        error('CDM_hybrid:Start','Parameter-name mapping requires unique labels. Supply an explicit Start for this layout change.');
    end
    [found,source] = ismember(O.ParameterNames,R.ParameterNames);
    b(found) = R.bhat(source(found));
else
    b = R.bhat(:);
end
end

function yes = sameFamily(R,O)
yes = true;
if ~isfield(R,'EstimOpt'), return; end
for name = {'NB','Truncated','LogSeries','LatentDistribution'}
    field = name{1};
    if isfield(R.EstimOpt,field)
        if strcmp(field,'LatentDistribution'), equal = strcmp(R.EstimOpt.(field),O.(field));
        else, equal = isequal(double(R.EstimOpt.(field)),double(O.(field)));
        end
        if ~equal, yes = false; return; end
    end
end
end

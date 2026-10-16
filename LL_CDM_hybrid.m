function [LL,score,Details] = LL_CDM_hybrid(b,I,O)
% LL_CDM_hybrid Per-row weighted JOINT LOG likelihood and analytic log score.
%   [LL,S,D] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt) uses resolved input/options
%   from CDM_hybrid (EvaluateOnly=true can prepare them without estimation).
%   Unlike LL_CDM, this function returns LOG likelihood, not its negative.
%   sum(LL) is maximized; S(i,:) = derivative of LL(i) wrt raw b.
%
%   Integration uses fixed quadrature over population-normal innovations.
%   Truncated=1 normalizes AFTER integrating and includes
%   selection's derivative. Missing indicators contribute a likelihood of one.
%   Population and selected EY integrate the prior latent distribution;
%   posterior scores are descriptive E[LV|Y,indicators,X], not regressors.
%   With LogSeries=true, the integrated denominator is an unnormalized
%   positive-selection weight, not a population participation probability.

if ~isfield(O,'HybridVersion') || O.HybridVersion~=1
    error('CDM_hybrid:Options','Use options resolved by CDM_hybrid (EvaluateOnly=true is available).');
end
b = b(:);
if numel(b)~=O.Params || ~isreal(b) || any(~isfinite(b))
    error('CDM_hybrid:Input','b must be a finite real vector matching EstimOpt.Params.');
end
n = numel(I.Y); L = O.NLatent; J = size(I.Xmea,2); P = O.Params;
logseries = isfield(O,'LogSeries') && O.LogSeries;
selected = isfield(O,'LatentDistribution') && strcmp(O.LatentDistribution,'selected');
Q = numel(O.QuadratureWeights); ix = O.Layout;
beta = b(ix.beta); gamma = reshape(b(ix.gamma),O.NVarStr,L);
delta = zeros(O.NVarLatent,L);
if O.LatentActive, delta = reshape(b(ix.delta),O.NVarLatent,L); end
loading = b(ix.loading); loading(O.Anchor) = exp(loading(O.Anchor));
cuts = cell(J,1); gaps = cell(J,1);
for j = 1:J
    raw = b(ix.thresholds{j});
    gaps{j} = [1;exp(raw(2:end))];
    cuts{j} = raw(1)+[0;cumsum(gaps{j}(2:end))];
end
mu = I.Xstr*gamma;
effects = I.Xlatent*delta;
z = cell(L,1);
eta = repmat(I.Xa*beta,1,Q);
for l = 1:L
    z{l} = mu(:,l)+O.QuadratureNodes(:,l)';
    eta = eta+z{l}.*effects(:,l);
end
logalpha = [];
if O.NB, logalpha = b(ix.logalpha); end
if logseries
    % Selected NB2 limit theta=1/alpha -> 0 at fixed v=log(alpha*mu):
    % f(y|z)/theta -> logistic(v)^y/y and Pr(Y>0|z)/theta -> softplus(v).
    % Cancel theta before integrating; the remaining denominator tilts the
    % population latent distribution by softplus(v), not by a probability.
    logp = -softplus(-eta);
    logCount = I.Y.*logp-log(I.Y);
    dc = I.Y.*logistic(-eta);
    da = zeros(size(eta));
    logSelect = log(softplus(eta));
    small = eta < -35;
    logSelect(small) = eta(small); % log(softplus(v)) = v to double precision
    dSelect = exp(logp-logSelect);
    daSelect = zeros(size(eta));
else
    [logCount,dc,da,logP0,dc0,da0] = countKernel(I.Y,eta,logalpha);
    logSelect = log(-expm1(logP0));
    ratio = exp(logP0-logSelect);
    dSelect = -weightedProduct(ratio,dc0);
    daSelect = -weightedProduct(ratio,da0);
end
logMeasurement = zeros(n,Q);
for j = 1:J
    lp = itemKernel(I.Xmea(:,j),loading(j)*z{O.ItemFactor(j)},cuts{j});
    logMeasurement = logMeasurement+lp;
end
if isfield(O,'QuadratureLogWeights')
    logw = O.QuadratureLogWeights(:)';
else
    logw = log(O.QuadratureWeights(:)');
end
logJoint = logCount+logMeasurement+logw;
if selected, logJoint = logJoint-logSelect; end
[joint,post] = logsum(logJoint);
[logSelection,selectionPost] = logsum(logSelect+logw);
LL = joint;
if O.Truncated && ~selected, LL = LL-logSelection; end
LL = I.W.*LL;

score = zeros(n,P);
if nargout>1 && all(isfinite(LL))
    % Subtract the denominator score averaged over its tilted population
    % distribution, which excludes the measurement indicators.
    sc = weightedProduct(post,dc);
    sa = weightedProduct(post,da);
    if selected
        % Here the normal law is already conditional on visiting. The
        % conditional count PMF's normalization is inside the same posterior.
        sc = sc-weightedProduct(post,dSelect);
        sa = sa-weightedProduct(post,daSelect);
    elseif O.Truncated
        sc = sc-weightedProduct(selectionPost,dSelect);
        sa = sa-weightedProduct(selectionPost,daSelect);
    end
    score(:,ix.beta) = sum(sc,2).*I.Xa;
    if O.NB, score(:,ix.logalpha) = sum(sa,2); end
    structural = zeros(n,L);
    for l = 1:L
        structural(:,l) = sum(sc,2).*effects(:,l);
        if O.LatentActive, score(:,ix.delta(:,l)) = sum(sc.*z{l},2).*I.Xlatent; end
    end
    for j = 1:J
        l = O.ItemFactor(j);
        [~,du,dt] = itemKernel(I.Xmea(:,j),loading(j)*z{l},cuts{j});
        weighted = post.*du;
        dl = sum(weighted.*z{l},2);
        if O.Anchor(j), dl = dl*loading(j); end
        score(:,ix.loading(j)) = dl;
        structural(:,l) = structural(:,l)+loading(j)*sum(weighted,2);
        for k = 1:numel(ix.thresholds{j})
            score(:,ix.thresholds{j}(k)) = sum(post.*dt{k},2)*gaps{j}(k);
        end
    end
    for l = 1:L, score(:,ix.gamma(:,l)) = structural(:,l).*I.Xstr; end
    score = score.*I.W;
end

if nargout>2
    [~,postItems] = logsum(logMeasurement+logw);
    [logSelectionItems,postItemsSelected] = logsum(logMeasurement+logSelect+logw);
    [logItems,~] = logsum(logMeasurement+logw);
    if selected, postItemsSelected = postItems; end
    Details = struct();
    Details.Parameters = struct('beta',beta,'gamma',gamma,'delta',delta, ...
        'loadings',loading,'thresholds',{cuts},'alpha',[]);
    Details.Parameters.LatentDistribution = 'population';
    if selected, Details.Parameters.LatentDistribution = 'selected'; end
    if O.NB, Details.Parameters.alpha = exp(logalpha); end
    Details.LatentMean = mu;
    Details.PosteriorMean = zeros(n,L);
    Details.PosteriorMeanIndicators = zeros(n,L);
    Details.PosteriorMeanIndicatorsSelected = zeros(n,L);
    for l = 1:L
        Details.PosteriorMean(:,l) = sum(post.*z{l},2);
        Details.PosteriorMeanIndicators(:,l) = sum(postItems.*z{l},2);
        Details.PosteriorMeanIndicatorsSelected(:,l) = sum(postItemsSelected.*z{l},2);
    end
    % E[exp(index)|X] is analytic under the normal factor distribution. In
    % the boundary family it equals lim E[Y|X]/theta, not a population mean.
    logExpectedExp = I.Xa*beta+sum(mu.*effects,2)+0.5*sum(effects.^2,2);
    [logExpItems,~] = logsum(logMeasurement+eta+logw);
    expGivenIndicators = exp(logExpItems-logItems);
    weightGivenIndicators = exp(logSelectionItems-logItems);
    Details.EYSelected = exp(logExpectedExp-logSelection);
    Details.EYGivenIndicatorsSelected = exp(logExpItems-logSelectionItems);
    Details.SelectionWeightIntegral = exp(logSelection);
    Details.SelectionWeightGivenIndicators = weightGivenIndicators;
    Details.LogSelectionWeightIntegral = logSelection;
    if logseries
        Details.EYPopulation = nan(n,1);
        Details.SelectionProbability = nan(n,1);
        Details.EYGivenIndicatorsPopulation = nan(n,1);
        Details.SelectionProbabilityGivenIndicators = nan(n,1);
        Details.LogSelectionProbability = nan(n,1);
        Details.Parameters.CountIndexMeaning = 'log(alpha*mu) in the selected NB2 boundary';
    else
        Details.EYPopulation = exp(logExpectedExp);
        Details.SelectionProbability = exp(logSelection);
        Details.EYGivenIndicatorsPopulation = expGivenIndicators;
        Details.SelectionProbabilityGivenIndicators = weightGivenIndicators;
        Details.LogSelectionProbability = logSelection;
        Details.Parameters.CountIndexMeaning = 'log(mu)';
    end
    Details.LogJointNumerator = joint;
    Details.LogIndicatorMarginalSelected = logSelectionItems-logSelection;
    Details.PredictionNote = ['EYPopulation=E[Y|X]; EYSelected=E[Y|X,Y>0]. ' ...
        'EYGivenIndicatorsPopulation conditions on items only in the population; ' ...
        'EYGivenIndicatorsSelected conditions on items and Y>0. PosteriorMean additionally conditions on observed Y.'];
    if logseries
        Details.PredictionNote = ['Selected NB2 log-series boundary: the count index is log(alpha*mu). ' ...
            'Only selected means are defined by this limiting model. Population means and participation probabilities are NaN. ' ...
            'SelectionWeightIntegral integrates softplus(index); it is not a probability. ' ...
            'PosteriorMean uses Y and items; PosteriorMeanIndicatorsSelected tilts the item posterior by softplus(index).'];
    end
    if selected
        [logSelectedMean,~] = logsum(eta-logSelect+logw);
        [logSelectedMeanItems,~] = logsum(logMeasurement+eta-logSelect+logw);
        Details.EYSelected = exp(logSelectedMean);
        Details.EYGivenIndicatorsSelected = exp(logSelectedMeanItems-logItems);
        Details.LogIndicatorMarginalSelected = logItems;
        Details.EYPopulation = nan(n,1);
        Details.SelectionProbability = nan(n,1);
        Details.EYGivenIndicatorsPopulation = nan(n,1);
        Details.SelectionProbabilityGivenIndicators = nan(n,1);
        Details.LogSelectionProbability = nan(n,1);
        Details.SelectionWeightIntegral = nan(n,1);
        Details.SelectionWeightGivenIndicators = nan(n,1);
        Details.LogSelectionWeightIntegral = nan(n,1);
        Details.PredictionNote = ['The normal latent distribution describes visitors. Conditional positive-count ' ...
            'PMFs are mixed inside the integral; posterior item-only scores therefore use the visitor-normal law. ' ...
            'EYSelected and EYGivenIndicatorsSelected are visitor means. Population means, participation probabilities ' ...
            'and population selection weights are not reported. LogIndicatorMarginalSelected allows conditional ' ...
            'count scoring by subtracting it from the unweighted joint log likelihood.'];
    end
end
end

function [lp,dc,da,lp0,dc0,da0] = countKernel(y,eta,logalpha)
% Log-domain vectorized Poisson/NB2 kernel. The existing LL_CDM uses a design
% matrix per evaluation; this form handles the entire latent quadrature grid
% and avoids NB2 gamma-function cancellation near the Poisson boundary.
lambda = exp(eta);
if isempty(logalpha)
    lp = y.*eta-gammaln(y+1)-lambda;
    dc = y-lambda;
    lp0 = -lambda;
    dc0 = -lambda;
    da = zeros(size(eta)); da0 = da;
else
    alpha = exp(logalpha); theta = exp(-logalpha);
    t = eta+logalpha;
    s = softplus(t);
    r = logistic(t);
    % gamma(theta+y)/gamma(theta) / theta^y = prod_{k=0}^{y-1}(1+alpha*k).
    % Counts are integers; cumulated log products keep small alpha accurate.
    % ponytail: cap the lookup at 100000; large counts use gamma/psi identities.
    if max(y)<=100000
        k = (0:max(y)-1)';
        coeff = [0;cumsum(log1p(alpha*k))];
        deriv = [0;cumsum((alpha*k)./(1+alpha*k))];
        c = coeff(y+1); a = deriv(y+1);
    else
        c = gammaln(theta+y)-gammaln(theta)-y*log(theta);
        a = y-theta*(psi(theta+y)-psi(theta));
    end
    sm = s-r;
    small = t < -9;
    u = exp(t(small));
    sm(small) = u.^2.*(0.5+u.*(-2/3+u.*(0.75-0.8*u)));
    lp = c-gammaln(y+1)+y.*eta-(y+theta).*s;
    dc = y-(y+theta).*r;
    da = a-y.*r+theta*sm;
    lp0 = -theta*s;
    dc0 = -theta*r;
    da0 = theta*sm;
end
end

function [lp,du,dt] = itemKernel(category,index,tau)
% Stable difference of logistic CDFs, including the two unbounded categories.
[n,Q] = size(index); K = numel(tau)+1;
observed = ~isnan(category);
upper = inf(n,1); lower = -inf(n,1);
iu = observed & category<K; il = observed & category>1;
upper(iu) = tau(category(iu)); lower(il) = tau(category(il)-1);
a = upper-index; b = lower-index;
lp = -softplus(-a)-softplus(b)+log(-expm1(lower-upper));
lp(~observed,:) = 0;
if nargout>1
    du = logistic(a)+logistic(b)-1;
    du(~observed,:) = 0;
    % d log(P)/d upper and lower boundary, before cumulative-gap transform.
    dup = exp(-softplus(-a)-softplus(a)-lp);
    dlo = -exp(-softplus(-b)-softplus(b)-lp);
    dup(~observed,:) = 0; dlo(~observed,:) = 0;
    dt = cell(K-1,1);
    dt{1} = -du; % common shift of every finite threshold
    for k = 2:K-1
        dt{k} = dup.*(iu & category>=k)+dlo.*(il & category-1>=k);
    end
end
end

function y = softplus(x)
y = max(x,0)+log1p(exp(-abs(x)));
end

function y = logistic(x)
y = exp(-softplus(-x));
end

function [value,weights] = logsum(x)
m = max(x,[],2);
v = exp(x-m);
total = sum(v,2);
value = m+log(total);
weights = v./total;
end

function out = weightedProduct(weight,value)
% Nodes with exactly zero numerical mass contribute zero even if a raw score
% overflows there. Nonfinite values at positive-mass nodes remain visible.
out = weight.*value;
out(weight==0) = 0;
end

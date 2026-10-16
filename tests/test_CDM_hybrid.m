function nFail = test_CDM_hybrid(verbose)
% Joint count/measurement likelihood checks, including population selection.
if nargin < 1, verbose = 1; end
cdm_test_paths;
rng(20260918);
nFail = 0;
checks = {@factorisation,@selection,@missingItems,@gradients, ...
    @multifactor,@nestedRestriction,@quadrature,@highOrderQuadrature, ...
    @normalGrid,@logseriesSelection,@logseriesGradient,@logseriesLimit, ...
    @invalidInputs,@syntheticEvidence,@estimation};
for c = 1:numel(checks)
    try
        detail = checks{c}();
        tag = 'ok';
    catch ME
        nFail = nFail + 1;
        tag = 'FAIL'; detail = ME.message;
    end
    if verbose, fprintf('  %-5s %-28s %s\n',tag,func2str(checks{c}),detail); end
end
if verbose, fprintf('test_CDM_hybrid: %d passed, %d failed\n',numel(checks)-nFail,nFail); end
end

function detail = factorisation
% With no latent count effect, the existing count model must be recovered.
[I,E] = example(24);
maxerr = 0;
for nb = 0:1
    for trunc = 0:1
        E.NB = nb; E.Truncated = trunc;
        I.Y = max(I.Y,trunc);
        R = CDM_hybrid(I,E);
        b = parameters(R,0);
        f = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
        C = struct('NB',nb,'NBtype',2,'Zinf',0,'Truncated',trunc, ...
            'Censored',0,'RealMin',0,'NVarA',2,'NVarNB',nb,'NVarZinf',0);
        bc = b(R.EstimOpt.Layout.beta);
        if nb, bc = [bc;b(R.EstimOpt.Layout.logalpha)]; end
        count = -LL_CDM(I.Y,I.Xa,ones(numel(I.Y),nb),zeros(numel(I.Y),0),C,bc);
        for r = 1:numel(I.Y)
            m = integral(@(z) measurement(z+latentMean(r,R,b),I.Xmea(r,:),b,R.EstimOpt).*normal(z),-10,10, ...
                'AbsTol',1e-12,'RelTol',1e-10);
            maxerr = max(maxerr,abs(f(r)-count(r)-log(m)));
        end
    end
end
assert(maxerr < 2e-7,'Zero latent count effect does not factorise: %.3g.',maxerr);
detail = sprintf('Poisson/NB2 x selected/unselected, max error %.2g',maxerr);
end

function detail = selection
% Summing over all counts AND all indicator patterns must give probability one.
% The selected indicator marginal must tilt with Pr(Y>0|latent), which catches
% the incorrect integral of a pointwise zero-truncated count likelihood.
patterns = dec2bin(0:7)-'0'+1;
nmax = 100;
I = struct('Y',repelem((1:nmax)',8),'Xa',ones(8*nmax,1), ...
    'Xmea',repmat(patterns,nmax,1));
E = struct('EvaluateOnly',true,'Display',0,'NNodes',64,'Truncated',1);
maxerr = 0;
for nb = 0:1
    E.NB = nb;
    R = CDM_hybrid(I,E);
    b = parameters(R,0.55);
    b(R.EstimOpt.Layout.beta) = -0.7;
    b(R.EstimOpt.Layout.thresholds{1}) = 0;
    f = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    mass = sum(exp(f));
    assert(abs(mass-1) < 2e-7,'Selected joint mass is %.12g.',mass);
    observed = sum(exp(f(I.Xmea(:,1)==2)));
    alpha = 0;
    if nb, alpha = exp(b(R.EstimOpt.Layout.logalpha)); end
    survival = @(z) positiveProbability(exp(-0.7+0.55*z),alpha);
    expected = integral(@(z) logistic(exp(b(R.EstimOpt.Layout.loading(1)))*z).*survival(z).*normal(z),-10,10, ...
        'AbsTol',1e-12,'RelTol',1e-10)/integral(@(z) survival(z).*normal(z),-10,10, ...
        'AbsTol',1e-12,'RelTol',1e-10);
    maxerr = max(maxerr,abs(observed-expected));
    assert(abs(observed-expected) < 2e-7,'Selected indicator marginal %.9g vs %.9g.',observed,expected);
    assert(observed > 0.54,'Selection failed to shift the indicator distribution.');
end
detail = sprintf('joint mass 1; independently integrated indicator marginal error %.2g',maxerr);
end

function detail = missingItems
[I,E] = example(24);
I.Xmea(1,2) = NaN;
I.Xmea(2,:) = NaN;
R = CDM_hybrid(I,E);
b = parameters(R,0.4);
f = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
maxerr = 0;
for r = 1:2
    alternatives = R.INPUT;
    missing = find(isnan(I.Xmea(r,:)));
    combos = cell(1,numel(missing));
    for j = 1:numel(missing), combos{j} = 1:R.EstimOpt.NCategories(missing(j)); end
    grids = cell(size(combos)); [grids{:}] = ndgrid(combos{:});
    mass = 0;
    for k = 1:numel(grids{1})
        for j = 1:numel(missing), alternatives.Xmea(r,missing(j)) = grids{j}(k); end
        fk = LL_CDM_hybrid(b,alternatives,R.EstimOpt);
        mass = mass + exp(fk(r));
    end
    maxerr = max(maxerr,abs(exp(f(r))-mass));
end
assert(maxerr < 1e-11,'Missing items do not equal marginalisation: %.3g.',maxerr);
detail = sprintf('one/all items missing retained and marginalised, error %.2g',maxerr);
end

function detail = gradients
[I,E] = example(28);
I.Xmea(2,2) = NaN; I.Xmea(5,:) = NaN;
I.W = linspace(0.5,1.5,numel(I.Y))';
I.Xlatent = [ones(numel(I.Y),1),mod((0:numel(I.Y)-1)',2)];
maxerr = 0;
for nb = 0:1
    for trunc = 0:1
        E.NB = nb; E.Truncated = trunc; I.Y = max(I.Y,trunc);
        R = CDM_hybrid(I,E); b = parameters(R,0.45);
        [f,g] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
        assert(all(isfinite(f)) && all(isfinite(g(:))),'Nonfinite likelihood/score.');
        numeric = zeros(size(g));
        for j = 1:numel(b)
            h = 1e-6*max(1,abs(b(j))); bp = b; bm = b;
            bp(j) = bp(j)+h; bm(j) = bm(j)-h;
            numeric(:,j) = (LL_CDM_hybrid(bp,R.INPUT,R.EstimOpt)-LL_CDM_hybrid(bm,R.INPUT,R.EstimOpt))/(2*h);
        end
        maxerr = max(maxerr,max(abs(g(:)-numeric(:))./max(1,abs(numeric(:)))));
        U = R.INPUT; U.W = ones(size(I.W));
        unweighted = LL_CDM_hybrid(b,U,R.EstimOpt);
        assert(max(abs(f-I.W.*unweighted)) < 1e-10,'Weights do not multiply respondent log likelihoods.');
    end
end
assert(maxerr < 2e-5,'Per-observation score relative error %.3g.',maxerr);
detail = sprintf('Poisson/NB2 x truncation, latent interactions + ordinal gaps + weights + NaNs, error %.2g',maxerr);
end

function detail = quadrature
[I,E] = example(24);
E.NB = 1; E.Truncated = 1; I.Y = max(I.Y,1); E.NNodes = 32;
R32 = CDM_hybrid(I,E); b = parameters(R32,0.5);
f32 = LL_CDM_hybrid(b,R32.INPUT,R32.EstimOpt);
E.NNodes = 64; R64 = CDM_hybrid(I,E);
f64 = LL_CDM_hybrid(b,R64.INPUT,R64.EstimOpt);
err = max(abs(f32-f64));
assert(err < 1e-6,'32 vs 64 nodes changed respondent likelihood by %.3g.',err);
detail = sprintf('32 vs 64 nodes: max respondent LL change %.2g',err);
end

function detail = multifactor
[I,E] = example(24);
I.Xmea = [I.Xmea,circshift(I.Xmea,3,1)];
E.MeaMatrix = blkdiag(ones(3,1),ones(3,1));
E.NNodes = 9; E.NB = 1; E.Truncated = 1; I.Y = max(I.Y,1);
R = CDM_hybrid(I,E); b = R.bhat; L = R.EstimOpt.Layout;
b(L.delta) = [0.3,-0.25]; b(L.gamma) = [0.2,-0.4];
b(L.loading) = [log(0.8);0.6;0.7;log(0.6);0.9;0.75];
[~,g] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
numeric = zeros(size(g));
for j = 1:numel(b)
    h = 1e-6*max(1,abs(b(j))); bp = b; bm = b;
    bp(j) = bp(j)+h; bm(j) = bm(j)-h;
    numeric(:,j) = (LL_CDM_hybrid(bp,R.INPUT,R.EstimOpt)-LL_CDM_hybrid(bm,R.INPUT,R.EstimOpt))/(2*h);
end
err = max(abs(g(:)-numeric(:))./max(1,abs(numeric(:))));
assert(err < 2e-5,'Two-factor score error %.3g.',err);
detail = sprintf('independent two-factor NB2 selected model: score error %.2g',err);
end

function detail = highOrderQuadrature
% High-order nodes need stable tail weights, beyond eigenvector precision.
[I,E] = example(24);
maxerr = 0;
for q = [256,512]
    E.NNodes = q; R = CDM_hybrid(I,E);
    z = R.EstimOpt.QuadratureNodes;
    w = R.EstimOpt.QuadratureWeights;
    err = max(abs([sum(w),sum(w.*z),sum(w.*z.^2),sum(w.*z.^4)]-[1,0,1,3]));
    maxerr = max(maxerr,err);
    assert(all(isfinite(R.EstimOpt.QuadratureLogWeights)),'High-order log weights lost tail nodes.');
    assert(err < 1e-9,'NNodes=%d normal moments inaccurate: %.3g.',q,err);
    assert(all(isfinite(R.LLdetailed)),'NNodes=%d likelihood is nonfinite.',q);
end
detail = sprintf('256/512 nodes reproduce normal moments 0..4; max error %.2g',maxerr);
end

function detail = normalGrid
% Sharp measurement likelihoods reproduce an independent adaptive integral.
[I,E] = example(24); I.Y = max(I.Y,1);
E.Truncated = 1; E.Quadrature = 'normal-grid'; E.GridLimit = 12;
gridError = 0; referenceError = 0;
for nb = 0:1
    E.NB = nb; E.NNodes = 128; A = CDM_hybrid(I,E);
    b = parameters(A,0.7); b(A.Layout.loading) = [log(6);6.3;4.8];
    fa = LL_CDM_hybrid(b,A.INPUT,A.EstimOpt);
    E.NNodes = 256; B = CDM_hybrid(I,E);
    fb = LL_CDM_hybrid(b,B.INPUT,B.EstimOpt);
    gridError = max(gridError,max(abs(fa-fb)));
    alpha = 0; if nb, alpha = exp(b(B.Layout.logalpha)); end
    for r = 1:4
        location = latentMean(r,B,b);
        mu = @(z) exp(I.Xa(r,:)*b(B.Layout.beta)+0.7*(z+location));
        numerator = integral(@(z) countProbability(I.Y(r),mu(z),alpha).* ...
            measurement(z+location,I.Xmea(r,:),b,B.EstimOpt).*normal(z),-12,12, ...
            'AbsTol',1e-13,'RelTol',1e-10);
        denominator = integral(@(z) positiveProbability(mu(z),alpha).*normal(z),-12,12, ...
            'AbsTol',1e-13,'RelTol',1e-10);
        referenceError = max(referenceError,abs(fb(r)-log(numerator/denominator)));
    end
end
assert(gridError < 1e-4,'Normal-grid 128/256 error %.3g.',gridError);
assert(referenceError < 1e-8,'Normal-grid vs adaptive integral error %.3g.',referenceError);
E.GridLimit = 3; T = CDM_hybrid(I,E);
mass = sum(T.EstimOpt.QuadratureWeights);
assert(abs(mass-erf(3/sqrt(2)))<2e-6 && mass<0.999, ...
    'Finite-grid normal mass was renormalized or integrated incorrectly.');
detail = sprintf('loadings ~6, Poisson/NB2: 128/256 %.2g; adaptive integral %.2g',gridError,referenceError);
end

function detail = logseriesSelection
patterns = dec2bin(0:7)-'0'+1; nmax = 150;
I = struct('Y',repelem((1:nmax)',8),'Xa',ones(8*nmax,1), ...
    'Xmea',repmat(patterns,nmax,1));
E = struct('EvaluateOnly',true,'Display',0,'NNodes',256,'Truncated',1, ...
    'LogSeries',true,'Quadrature','normal-grid');
R = CDM_hybrid(I,E); b = parameters(R,0.55);
b(R.Layout.beta) = -0.7; b(R.Layout.thresholds{1}) = 0;
f = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
mass = sum(exp(f));
observed = sum(exp(f(I.Xmea(:,1)==2)));
tilt = @(z) log1p(exp(-0.7+0.55*z));
denominator = integral(@(z) tilt(z).*normal(z),-12,12,'AbsTol',1e-12);
expected = integral(@(z) logistic(0.9*z).*tilt(z).*normal(z),-12,12, ...
    'AbsTol',1e-12)/denominator;
assert(abs(mass-1) < 1e-8,'Log-series joint mass %.12g.',mass);
assert(abs(observed-expected)<1e-8,'Log-series selected indicator distribution is incorrect.');
assert(observed>0.54,'Log-series indicator selection tilt was omitted.');
detail = sprintf('joint mass %.12g; independent tilted indicator error %.2g',mass,abs(observed-expected));
end

function detail = logseriesGradient
[I,E] = example(28); I.Y = max(I.Y,1);
I.Xmea(2,2) = NaN; I.Xmea(5,:) = NaN;
I.W = linspace(0.5,1.5,numel(I.Y))';
I.Xlatent = [ones(numel(I.Y),1),mod((0:numel(I.Y)-1)',2)];
E.LogSeries = true; E.Truncated = 1; E.Quadrature = 'normal-grid'; E.NNodes = 128;
R = CDM_hybrid(I,E); b = parameters(R,0.45);
[f,g] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
assert(all(isfinite(f)) && all(isfinite(g(:))),'Log-series nonfinite likelihood/score.');
numeric = zeros(size(g));
for j = 1:numel(b)
    h = 1e-6*max(1,abs(b(j))); bp = b; bm = b;
    bp(j) = bp(j)+h; bm(j) = bm(j)-h;
    numeric(:,j) = (LL_CDM_hybrid(bp,R.INPUT,R.EstimOpt)-LL_CDM_hybrid(bm,R.INPUT,R.EstimOpt))/(2*h);
end
err = max(abs(g(:)-numeric(:))./max(1,abs(numeric(:))));
assert(err<2e-5,'Log-series gradient error %.3g.',err);
detail = sprintf('structural, threshold gaps, interactions, weights, NaNs: score error %.2g',err);
end

function detail = logseriesLimit
[I,E] = example(24); I.Y = max(I.Y,1);
E.LogSeries = true; E.Truncated = 1; E.Quadrature = 'normal-grid'; E.NNodes = 256;
R = CDM_hybrid(I,E); b = parameters(R,0.55);
[f,~,D] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
assert(all(isnan(D.EYPopulation)) && all(isnan(D.SelectionProbability)) && ...
    all(isnan(D.EYGivenIndicatorsPopulation)) && all(isnan(D.SelectionProbabilityGivenIndicators)), ...
    'Boundary model reported identified population predictions.');
assert(all(isfinite(D.EYSelected)) && all(isfinite(D.EYGivenIndicatorsSelected)), ...
    'Boundary model omitted selected predictions.');
errs = zeros(2,1); predictionError = 0;
for j = 1:2
    a = [12,18]; a = a(j);
    F = E; F.LogSeries = false; F.NB = 1;
    N = CDM_hybrid(I,F); bn = parameters(N,0.55);
    bn(N.Layout.logalpha) = a;
    bn(N.Layout.beta(1)) = bn(N.Layout.beta(1))-a;
    [fn,~,Dn] = LL_CDM_hybrid(bn,N.INPUT,N.EstimOpt);
    errs(j) = max(abs(fn-f));
    if j==2
        predictionError = max([abs(Dn.EYSelected-D.EYSelected)./D.EYSelected; ...
            abs(Dn.EYGivenIndicatorsSelected-D.EYGivenIndicatorsSelected)./D.EYGivenIndicatorsSelected]);
    end
end
assert(errs(2)<1e-6 && errs(2)<errs(1)/100,'Selected NB did not converge to log-series: %s.',mat2str(errs));
assert(predictionError<1e-6,'Selected NB/log-series prediction relative difference %.3g.',predictionError);
detail = sprintf('NB logalpha12/18 LL errors %.2g/%.2g; selected EY error %.2g',errs(1),errs(2),predictionError);
end

function p = countProbability(y,mu,alpha)
if alpha==0
    p = exp(y*log(mu)-mu-gammaln(y+1));
else
    theta = 1/alpha;
    p = exp(gammaln(y+theta)-gammaln(theta)-gammaln(y+1)+ ...
        y*log(mu/theta)-(y+theta)*log1p(mu/theta));
end
end

function detail = nestedRestriction
[I,E] = example(24);
E.NB = 1; E.Truncated = 1; I.Y = max(I.Y,1);
R = CDM_hybrid(I,E); b = parameters(R,0);
E.LatentActive = false;
N = CDM_hybrid(I,E);
[found,ix] = ismember(N.ParameterNames,R.ParameterNames);
assert(all(found),'Nested model changed shared parameter names.');
assert(numel(b)-numel(N.bhat)==numel(R.Layout.delta),'Nested model removed the wrong parameter count.');
[f,g] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
[fn,gn] = LL_CDM_hybrid(b(ix),N.INPUT,N.EstimOpt);
err = max([abs(f-fn);max(abs(g(:,ix)-gn),[],2)]);
assert(err < 1e-10,'LatentActive=false does not nest zero delta: %.3g.',err);
detail = sprintf('zero delta and LatentActive=false identical, error %.2g',err);
end

function detail = invalidInputs
[I,E] = example(24);
F = E; F.NCT = 2; rejected(I,F,{'CDM_hybrid:Unsupported'});
F = E; F.Truncated = 2; rejected(I,F,{'CDM_hybrid:Unsupported','CDM_hybrid:Options'});
J = I; J.Xstr = ones(size(I.Y)); rejected(J,E,{'CDM_hybrid:Identification'});
J = I; J.Xstr = [double(I.Xa(:,2)<0),double(I.Xa(:,2)>=0)]; rejected(J,E,{'CDM_hybrid:Identification'});
J = I; J.Xa(:,3) = J.Xa(:,2); rejected(J,E,{'CDM:RankDeficientDesign','CDM_hybrid:Identification'});
J = I; J.Xmea(1,1) = 0; rejected(J,E,{'CDM_hybrid:Input','CDM_hybrid:Identification'});
J = I; J.Xmea = J.Xmea(:,1:2); rejected(J,E,{'CDM_hybrid:Identification'});
J = I; J.Xmea(:,1) = NaN; rejected(J,E,{'CDM_hybrid:Identification','CDM_hybrid:Input'});
F = E; F.LogSeries = true; rejected(I,F,{'CDM_hybrid:Options'});
F.Truncated = 1; F.NB = 1; rejected(I,F,{'CDM_hybrid:Options'});
F = E; F.Quadrature = 'unknown'; rejected(I,F,{'CDM_hybrid:Options'});
F = E; F.GridLimit = 0; rejected(I,F,{'CDM_hybrid:Options'});
detail = 'panel, stratification, rank, indicators, log-series and quadrature options';
end

function detail = syntheticEvidence
% A correctly oriented shared factor must improve the joint likelihood for
% independent Poisson AND overdispersed synthetic samples at known parameters.
gains = zeros(2,1);
for nb = 0:1
    [I,E] = synthetic(1800,nb,20260918+nb);
    E.EvaluateOnly = true; E.NNodes = 40;
    R = CDM_hybrid(I,E); b = truth(R);
    f = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    null = b; null(R.EstimOpt.Layout.delta) = 0;
    opposite = b; opposite(R.EstimOpt.Layout.delta) = -0.55;
    gains(nb+1) = sum(f)-sum(LL_CDM_hybrid(null,R.INPUT,R.EstimOpt));
    assert(gains(nb+1)>30,'Synthetic NB=%d did not detect shared factor, gain %.3g.',nb,gains(nb+1));
    assert(sum(f)>sum(LL_CDM_hybrid(opposite,R.INPUT,R.EstimOpt))+50,'Wrong factor sign fits too well.');
end
detail = sprintf('Poisson/NB2 LL gains over zero factor: %.1f / %.1f',gains);
end

function detail = estimation
[I,E] = synthetic(1800,0,782311);
E.NNodes = 24; E.EvaluateOnly = true;
P = CDM_hybrid(I,E); target = truth(P);
E.Start = target; E.Start(P.EstimOpt.Layout.delta) = 0.15;
E.EvaluateOnly = false;
O = optimoptions('fminunc','Display','off','Algorithm','quasi-newton', ...
    'MaxIterations',300,'MaxFunctionEvaluations',10000,'OptimalityTolerance',1e-5);
R = CDM_hybrid(I,E,O);
L = R.EstimOpt.Layout;
assert(all(isfinite(R.bhat)) && R.exitf>0,'Hybrid estimation did not converge.');
assert(max(abs(R.bhat(L.beta)-target(L.beta))) < 0.2,'Count coefficient recovery failed.');
assert(abs(R.bhat(L.delta)-0.55)<0.2,'Latent count coefficient recovery failed.');
assert(abs(R.bhat(L.gamma)-0.35)<0.3,'Structural coefficient recovery failed.');
assert(abs(R.LL-sum(R.LLdetailed)) < 1e-8,'Detailed likelihood does not add up.');
assert(R.LL >= sum(LL_CDM_hybrid(target,R.INPUT,R.EstimOpt))-1e-5,'Fit worse than known generating parameters.');
detail = sprintf('Poisson fit: beta %s; latent effect %.3f (true .55)',mat2str(R.bhat(L.beta)',3),R.bhat(L.delta));
end

function [I,E] = example(n)
x = linspace(-1,1,n)';
I = struct('Y',mod((0:n-1)',5),'Xa',[ones(n,1),x], ...
    'Xstr',sin(2*x),'Xlatent',ones(n,1), ...
    'Xmea',[mod((0:n-1)',2)+1,mod((0:n-1)',3)+1,mod((0:n-1)',4)+1]);
E = struct('EvaluateOnly',true,'Display',0,'NNodes',64);
end

function b = parameters(R,delta)
L = R.EstimOpt.Layout; b = R.bhat;
b(L.beta) = [0.3;repmat(-0.25,numel(L.beta)-1,1)];
b(L.logalpha) = log(0.35);
b(L.delta) = delta; b(L.gamma) = 0.3;
b(L.loading) = [log(0.9);0.7;1.1];
for j = 1:numel(L.thresholds)
    ix = L.thresholds{j}; b(ix) = [-0.4;repmat(log(1.1),numel(ix)-1,1)];
end
end

function p = measurement(z,items,b,E)
L = E.Layout; p = ones(size(z));
for j = 1:numel(items)
    if isnan(items(j)), continue; end
    loading = b(L.loading(j));
    if E.Anchor(j), loading = exp(loading); end
    raw = b(L.thresholds{j}); thresholds = [raw(1);raw(1)+cumsum(exp(raw(2:end)))];
    k = items(j); upper = ones(size(z)); lower = zeros(size(z));
    if k <= numel(thresholds), upper = logistic(thresholds(k)-loading*z); end
    if k > 1, lower = logistic(thresholds(k-1)-loading*z); end
    p = p.*(upper-lower);
end
end

function m = latentMean(r,R,b)
m = R.INPUT.Xstr(r,:)*b(R.EstimOpt.Layout.gamma);
end

function p = positiveProbability(mu,alpha)
if alpha == 0, p = -expm1(-mu);
else, p = -expm1(-log1p(alpha*mu)/alpha); end
end

function p = logistic(x)
p = 1./(1+exp(-x));
end

function p = normal(x)
p = exp(-x.^2/2)/sqrt(2*pi);
end

function rejected(I,E,ids)
try
    CDM_hybrid(I,E);
catch ME
    assert(any(strcmp(ME.identifier,ids)),'Wrong error for invalid input: %s (%s).',ME.identifier,ME.message);
    return
end
error('test_CDM_hybrid:AcceptedInvalid','Invalid hybrid specification was accepted.');
end

function [I,E] = synthetic(n,nb,seed)
rng(seed);
x = 2*rand(n,1)-1; s = 2*rand(n,1)-1;
z = 0.35*s+randn(n,1);
mu = exp(0.45-0.35*x+0.55*z);
if nb, y = nbinrnd(1/0.35,(1/0.35)./(1/0.35+mu));
else, y = poissrnd(mu); end
loads = [1,0.9,1.2]; thresholds = [-0.4,0.1,0.5];
items = 1+double(rand(n,3)>logistic(thresholds-z*loads));
I = struct('Y',y,'Xa',[ones(n,1),x],'Xstr',s,'Xmea',items);
E = struct('NB',nb,'Display',0,'NNodes',32);
end

function b = truth(R)
L = R.EstimOpt.Layout; b = R.bhat;
b(L.beta) = [0.45;-0.35]; b(L.logalpha) = log(0.35);
b(L.delta) = 0.55; b(L.gamma) = 0.35;
b(L.loading) = [0;0.9;1.2];
thresholds = [-0.4,0.1,0.5];
for j = 1:3, b(L.thresholds{j}) = thresholds(j); end
end

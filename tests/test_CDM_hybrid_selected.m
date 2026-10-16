function nFail = test_CDM_hybrid_selected(verbose)
% Visitor-normal sensitivity: conditional normalization, scores and means.
if nargin<1, verbose = 1; end
cdm_test_paths;
nFail = 0;
checks = {@normalization,@gradients,@zeroLink,@predictions,@options};
for k = 1:numel(checks)
    try
        detail = checks{k}(); tag = 'ok';
    catch ME
        nFail = nFail+1; detail = ME.message; tag = 'FAIL';
    end
    if verbose, fprintf('  %-5s %-28s %s\n',tag,func2str(checks{k}),detail); end
end
if verbose, fprintf('test_CDM_hybrid_selected: %d passed, %d failed\n',numel(checks)-nFail,nFail); end
end

function detail = normalization
patterns = dec2bin(0:7)-'0'+1;
I = struct('Y',repelem((1:120)',8),'Xa',ones(960,1),'Xmea',repmat(patterns,120,1));
maxerr = 0;
for family = 0:2
    O = settings(family); R = CDM_hybrid(I,O); b = values(R);
    [f,~,D] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    assert(abs(sum(exp(f))-1)<2e-8,'Selected-normal joint mass is %.12g.',sum(exp(f)));
    for r = [1,8,51,128]
        independent = integral(@(z) countprob(I.Y(r),exp(b(R.Layout.beta)+0.35*z),family).* ...
            measurement(z,I.Xmea(r,:),b,R.EstimOpt).*normal(z),-12,12, ...
            'AbsTol',1e-13,'RelTol',1e-10);
        maxerr = max(maxerr,abs(exp(f(r))-independent));
    end
    % Summing over counts removes the conditional count model completely.
    for j = 1:8
        itemMass = sum(exp(f(j:8:end)));
        independent = integral(@(z) measurement(z,patterns(j,:),b,R.EstimOpt).*normal(z),-12,12);
        assert(abs(itemMass-independent)<2e-8,'Visitor-normal indicator marginal is incorrectly selection-tilted.');
        assert(abs(exp(D.LogIndicatorMarginalSelected(j))-independent)<1e-9,'Incorrect item marginal output.');
    end
end
assert(maxerr<1e-9,'Adaptive integral mismatch %.3g.',maxerr);
detail = sprintf('Poisson/NB2/log-series, joint mass and item marginal correct, integral error %.2g',maxerr);
end

function detail = gradients
[I,~] = example;
I.Xmea(2,2) = NaN; I.Xmea(5,:) = NaN;
I.W = linspace(0.5,1.5,numel(I.Y))';
I.Xlatent = [ones(numel(I.Y),1),mod((0:numel(I.Y)-1)',2)];
maxerr = 0;
for family = 0:2
    R = CDM_hybrid(I,settings(family)); b = values(R);
    [~,g] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    numeric = zeros(size(g));
    for j = 1:numel(b)
        step = 1e-6*max(1,abs(b(j))); bp=b; bm=b;
        bp(j)=bp(j)+step; bm(j)=bm(j)-step;
        numeric(:,j) = (LL_CDM_hybrid(bp,R.INPUT,R.EstimOpt)-LL_CDM_hybrid(bm,R.INPUT,R.EstimOpt))/(2*step);
    end
    maxerr = max(maxerr,max(abs(g(:)-numeric(:))./max(1,abs(numeric(:)))));
end
assert(maxerr<2e-5,'Selected-normal score relative error %.3g.',maxerr);
detail = sprintf('all three families, interactions/structural/items/weights/missing; error %.2g',maxerr);
end

function detail = zeroLink
[I,~] = example;
maxerr = 0;
for family = 0:2
    O = settings(family); R = CDM_hybrid(I,O); b = values(R); b(R.Layout.delta)=0;
    fSelected = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    O = R.EstimOpt; O.LatentDistribution = 'population';
    fPopulation = LL_CDM_hybrid(b,R.INPUT,O);
    maxerr = max(maxerr,max(abs(fSelected-fPopulation)));
end
assert(maxerr<1e-10,'Zero-link models differ by %.3g.',maxerr);
detail = sprintf('population and visitor normal coincide at delta=0, error %.2g',maxerr);
end

function detail = predictions
[I,~] = example;
maxerr = 0;
for family = 0:2
    R = CDM_hybrid(I,settings(family)); b = values(R);
    [~,~,D] = LL_CDM_hybrid(b,R.INPUT,R.EstimOpt);
    assert(all(isnan(D.EYPopulation)) && all(isnan(D.SelectionProbability)), ...
        'Selected-only model reports a population estimand.');
    assert(max(abs(D.PosteriorMeanIndicatorsSelected-D.PosteriorMeanIndicators),[],'all')<1e-12, ...
        'Selected item posterior has an extra selection tilt.');
    for r = [1,8,17]
        meanLV = I.Xstr(r,:)*b(R.Layout.gamma);
        index = I.Xa(r,:)*b(R.Layout.beta);
        h = b(R.Layout.delta);
        avg = integral(@(e) conditionalMean(exp(index+h*(meanLV+e)),family).*normal(e), ...
            -12,12,'AbsTol',1e-11,'RelTol',1e-9);
        denominator = integral(@(e) measurement(meanLV+e,I.Xmea(r,:),b,R.EstimOpt).*normal(e),-12,12);
        conditional = integral(@(e) conditionalMean(exp(index+h*(meanLV+e)),family).* ...
            measurement(meanLV+e,I.Xmea(r,:),b,R.EstimOpt).*normal(e),-12,12)/denominator;
        maxerr = max([maxerr,abs(avg-D.EYSelected(r)),abs(conditional-D.EYGivenIndicatorsSelected(r))]);
    end
end
assert(maxerr<1e-7,'Selected mean mismatch %.3g.',maxerr);
detail = sprintf('prior/items-conditioned selected means match independent integrals, error %.2g',maxerr);
end

function detail = options
[I,O] = example;
for invalid = {struct('LatentDistribution','unknown'),struct('LatentDistribution','selected','Truncated',0)}
    bad=O; fields=fieldnames(invalid{1});
    for j=1:numel(fields), bad.(fields{j})=invalid{1}.(fields{j}); end
    rejected=false;
    try, CDM_hybrid(I,bad); catch ME, rejected=strcmp(ME.identifier,'CDM_hybrid:Options'); end
    assert(rejected,'Invalid latent-distribution option was not rejected.');
end
R=CDM_hybrid(I,O); legacy=rmfield(R.EstimOpt,'LatentDistribution');
assert(isequal(LL_CDM_hybrid(R.bhat,R.INPUT,legacy),R.LLdetailed),'Legacy population options changed.');
O.LatentDistribution='selected'; R=CDM_hybrid(I,O);
assert(endsWith(R.Model,'_SELECTED_LATENT'),'Model does not distinguish latent sampling assumption.');
detail='invalid modes rejected, old saved options supported, distinct model label';
end

function O = settings(family)
O = struct('EvaluateOnly',true,'Display',0,'NNodes',128,'Quadrature','normal-grid', ...
    'GridLimit',12,'Truncated',1,'LatentDistribution','selected','NB',family==1,'LogSeries',family==2);
end

function [I,O] = example
n=32; x=linspace(-1,1,n)';
I = struct('Y',mod((0:n-1)',5)+1,'Xa',[ones(n,1),x],'Xstr',sin(2*x), ...
    'Xmea',[mod((0:n-1)',2)+1,mod((0:n-1)',3)+1,mod((0:n-1)',4)+1]);
O=settings(0); O.LatentDistribution='population';
end

function b = values(R)
ix=R.Layout; b=R.bhat;
b(ix.beta)=[0.2;repmat(-0.3,numel(ix.beta)-1,1)];
b(ix.delta)=0.35; b(ix.gamma)=0.2; b(ix.logalpha)=log(0.35);
b(ix.loading)=[log(0.9);0.7;1.1];
for j=1:numel(ix.thresholds)
    k=ix.thresholds{j}; b(k)=[-0.4;repmat(log(1.1),numel(k)-1,1)];
end
end

function p = measurement(z,item,b,O)
p=ones(size(z));
for j=1:numel(item)
    if isnan(item(j)), continue; end
    loading=b(O.Layout.loading(j)); if O.Anchor(j), loading=exp(loading); end
    raw=b(O.Layout.thresholds{j}); tau=[raw(1);raw(1)+cumsum(exp(raw(2:end)))];
    upper=ones(size(z)); lower=zeros(size(z));
    if item(j)<=numel(tau), upper=1./(1+exp(-(tau(item(j))-loading*z))); end
    if item(j)>1, lower=1./(1+exp(-(tau(item(j)-1)-loading*z))); end
    p=p.*(upper-lower);
end
end

function p = countprob(y,mu,family)
if family==2
    p=(mu./(1+mu)).^y./(y*log1p(mu));
elseif family==1
    theta=1/0.35;
    p=exp(gammaln(y+theta)-gammaln(theta)-gammaln(y+1)+theta*log(theta)+ ...
        y*log(mu)-(y+theta)*log(theta+mu))./(-expm1(-theta*log1p(mu/theta)));
else
    p=exp(y*log(mu)-mu-gammaln(y+1))./(-expm1(-mu));
end
end

function out = conditionalMean(mu,family)
if family==2, out=mu./log1p(mu);
elseif family==1, out=mu./(-expm1(-log1p(0.35*mu)/0.35));
else, out=mu./(-expm1(-mu)); end
end

function p=normal(z)
p=exp(-z.^2/2)/sqrt(2*pi);
end

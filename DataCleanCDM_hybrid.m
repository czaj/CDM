function [I,Results,O,OptimOpt] = DataCleanCDM_hybrid(I,O)
% DataCleanCDM_hybrid validates hybrid count data and resolves model options.
%   [INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM_hybrid(INPUT,EstimOpt)
%   follows the four-output CDM/DCE preparation convention. Results is an
%   empty accumulator: no optimized null-model likelihood is fabricated.
%   EstimOpt.InitialValues contains the identified raw starting vector.
%
%   Count/design rows must already be complete; NaN ordinal indicators are
%   retained and integrated out. Supplied weights are preserved exactly.
%   Include a constant in Xa but never in Xstr. See CDM_hybrid for options.
%
%   SaveXLS and SaveTxtOutput default to Display~=0. EvaluateOnly suppresses
%   automatic output setup. Export uses the companion Tools repository;
%   quiet numerical estimation has no Tools dependency.
%
%   See also CDM_hybrid, LL_CDM_hybrid, LL_CDM_hybrid_MATlike.
if nargin<2 || isempty(O), O=struct(); end
if ~isstruct(I) || ~isstruct(O)
    error('CDM_hybrid:Input','INPUT and EstimOpt must be structures.');
end
defaults = struct('NB',0,'NBtype',2,'Truncated',0,'Censored',0,'Zinf',0, ...
    'NCT',1,'LatentActive',true,'NNodes',15,'MaxQuadratureNodes',4096, ...
    'Display',1,'EvaluateOnly',false,'ComputeHessian',true,'RobustStd',false, ...
    'LogSeries',false,'Quadrature','hermite','GridLimit',12,'LatentDistribution','population');
fields = fieldnames(defaults);
for j = 1:numel(fields)
    if ~isfield(O,fields{j}) || isempty(O.(fields{j})), O.(fields{j}) = defaults.(fields{j}); end
end
for name = {'NB','Truncated','LatentActive','Display','EvaluateOnly','ComputeHessian','RobustStd','LogSeries'}
    value = O.(name{1});
    if ~isscalar(value) || ~isreal(value) || ~ismember(value,[0,1])
        error('CDM_hybrid:Options','%s must be 0 or 1.',name{1});
    end
end
if O.LogSeries && (O.NB~=0 || O.Truncated~=1)
    error('CDM_hybrid:Options','LogSeries=true requires NB=0 and Truncated=1.');
end
if ~(ischar(O.LatentDistribution) || (isstring(O.LatentDistribution) && isscalar(O.LatentDistribution))) || ...
        ~ismember(string(O.LatentDistribution),["population","selected"])
    error('CDM_hybrid:Options','LatentDistribution must be ''population'' or ''selected''.');
end
O.LatentDistribution = char(O.LatentDistribution);
if strcmp(O.LatentDistribution,'selected') && O.Truncated~=1
    error('CDM_hybrid:Options','LatentDistribution=''selected'' requires Truncated=1.');
end
if ~(ischar(O.Quadrature) || (isstring(O.Quadrature) && isscalar(O.Quadrature))) || ...
        ~ismember(string(O.Quadrature),["hermite","normal-grid"])
    error('CDM_hybrid:Options','Quadrature must be ''hermite'' or ''normal-grid''.');
end
O.Quadrature = char(O.Quadrature);
if ~isscalar(O.GridLimit) || ~isreal(O.GridLimit) || ~isfinite(O.GridLimit) || O.GridLimit<=0
    error('CDM_hybrid:Options','GridLimit must be a positive finite scalar.');
end
if O.NBtype ~= 2 || O.NCT ~= 1 || O.Censored ~= 0 || O.Zinf ~= 0
    error('CDM_hybrid:Unsupported','Only cross-sectional, uncensored Poisson/NB2 without zero inflation is supported.');
end
for name = {'ConstVarActive','Hurdle','probit','NumGrad'}
    if isfield(O,name{1}) && any(O.(name{1})(:) ~= 0)
        error('CDM_hybrid:Unsupported','The option %s is unsupported; remove it.',name{1});
    end
end
if isfield(O,'BActive') && any(O.BActive(:) ~= 1)
    error('CDM_hybrid:Unsupported','Use LatentActive=false for the nested count restriction; BActive restrictions are unsupported.');
end
if ~isfield(I,'Y') || ~isnumeric(I.Y) || ~isreal(I.Y) || ~isvector(I.Y) || isempty(I.Y)
    error('CDM_hybrid:Input','Y must be a nonempty numeric count vector.');
end
I.Y = double(I.Y(:)); n = numel(I.Y);
if any(~isfinite(I.Y) | I.Y<0 | I.Y~=fix(I.Y)) || (O.Truncated && any(I.Y==0))
    error('CDM_hybrid:Input','Counts must be finite nonnegative integers; Truncated=1 requires positive counts.');
end
if ~isfield(I,'Xa') || isempty(I.Xa), error('CDM_hybrid:Input','Provide count covariates Xa with your own intercept.'); end
if ~isfield(I,'Xmea') || isempty(I.Xmea), error('CDM_hybrid:Input','Provide ordinal indicators Xmea.'); end
if ~isfield(I,'Xstr') || isempty(I.Xstr), I.Xstr = zeros(n,0); end
if ~isfield(I,'Xlatent') || isempty(I.Xlatent), I.Xlatent = ones(n,1); end
if ~isfield(I,'W') || isempty(I.W), I.W = ones(n,1); end
I.W = I.W(:);
for name = {'Xa','Xstr','Xlatent','Xmea','W'}
    value = I.(name{1});
    if ~isnumeric(value) || ~isreal(value) || ~ismatrix(value) || size(value,1)~=n
        error('CDM_hybrid:Input','%s must be a real numeric matrix with n rows.',name{1});
    end
    if strcmp(name{1},'Xmea'), bad = isinf(value); else, bad = ~isfinite(value); end
    if any(bad(:)), error('CDM_hybrid:Input','Invalid nonfinite values in %s.',name{1}); end
    I.(name{1}) = double(value);
end
if size(I.W,2)~=1 || any(I.W<=0), error('CDM_hybrid:Input','W must contain one positive finite weight per row.'); end
if rank(I.Xa)<size(I.Xa,2) || rank(I.Xstr)<size(I.Xstr,2) || rank(I.Xlatent)<size(I.Xlatent,2)
    error('CDM_hybrid:Identification','Xa, Xstr and Xlatent must each have full column rank.');
end
if ~isempty(I.Xstr) && rank([ones(n,1),I.Xstr])==rank(I.Xstr)
    error('CDM_hybrid:Identification','Xstr cannot contain a constant: latent location is fixed at zero for Xstr=0.');
end
for name = {'Xnb','Xzinf'}
    if isfield(I,name{1}) && ~isempty(I.(name{1}))
        error('CDM_hybrid:Unsupported','%s is unsupported; NB2 has a single constant log dispersion.',name{1});
    end
end
if isfield(I,'MissingInd') && any(I.MissingInd(:)~=0)
    error('CDM_hybrid:Unsupported','Remove missing count/design rows before fitting; NaN measurement items are supported.');
end
J = size(I.Xmea,2);
if ~isfield(O,'MeaMatrix') || isempty(O.MeaMatrix), O.MeaMatrix = ones(J,1); end
M = O.MeaMatrix;
if ~(isnumeric(M)||islogical(M)) || ~isreal(M) || size(M,1)~=J || ...
        any(~ismember(M(:),[0,1])) || any(sum(M,2)~=1) || any(sum(M,1)<3)
    error('CDM_hybrid:Identification','MeaMatrix must assign each item to exactly one factor, with at least three items per factor.');
end
L = size(M,2);
if L<1 || L>3, error('CDM_hybrid:Unsupported','One to three independent latent factors are supported.'); end
if isfield(O,'NLatent') && ~isempty(O.NLatent) && O.NLatent~=L
    error('CDM_hybrid:Options','NLatent must match MeaMatrix columns.');
end
[~,factor] = max(M,[],2);
anchor = false(J,1);
for l = 1:L, anchor(find(factor==l,1)) = true; end
K = zeros(J,1);
for j = 1:J
    categories = unique(I.Xmea(~isnan(I.Xmea(:,j)),j));
    if numel(categories)<2 || any(categories' ~= 1:numel(categories))
        error('CDM_hybrid:Identification','Item %d must have at least two observed consecutive integer categories starting at 1.',j);
    end
    K(j) = numel(categories);
end
if isfield(O,'NCategories') && ~isempty(O.NCategories) && ~isequal(O.NCategories(:),K)
    error('CDM_hybrid:Identification','NCategories must match the observed categories; empty categories imply unidentified threshold gaps.');
end
nodeCap = 512;
if strcmp(O.Quadrature,'normal-grid'), nodeCap = O.MaxQuadratureNodes; end
if ~isscalar(O.NNodes) || O.NNodes<3 || O.NNodes>nodeCap || O.NNodes~=fix(O.NNodes) || ...
        ~isscalar(O.MaxQuadratureNodes) || ~isfinite(O.MaxQuadratureNodes) || O.MaxQuadratureNodes<1 || O.NNodes^L>O.MaxQuadratureNodes
    error('CDM_hybrid:Options','NNodes must be an integer >=3 (Hermite <=512); NNodes^NLatent cannot exceed MaxQuadratureNodes.');
end
if strcmp(O.Quadrature,'normal-grid')
    spacing = 2*O.GridLimit/O.NNodes;
    nodes = -O.GridLimit+spacing*((1:O.NNodes)'-0.5);
    logweights = log(spacing)-0.5*nodes.^2-0.5*log(2*pi);
else
    % Golub-Welsch for the standard-normal measure: no quadrature dependency.
    T = diag(sqrt(1:O.NNodes-1),1)+diag(sqrt(1:O.NNodes-1),-1);
    nodes = sort(eig(T));
    % At high orders tiny first eigenvector entries lose relative accuracy.
    % w_i = 1/(N*p_{N-1}(x_i)^2), with orthonormal Hermite p, gives stable log
    % weights. Rescaling the recurrence prevents squaring/recurrence overflow.
    previous = ones(size(nodes)); current = nodes; logscale = zeros(size(nodes));
    for degree = 1:O.NNodes-2
        next = (nodes.*current-sqrt(degree)*previous)/sqrt(degree+1);
        previous = current; current = next;
        large = max(abs(previous),abs(current))>1e100;
        previous(large) = previous(large)/1e100;
        current(large) = current(large)/1e100;
        logscale(large) = logscale(large)+log(1e100);
    end
    logweights = -log(O.NNodes)-2*(log(abs(current))+logscale);
    logweights = logweights-log(sum(exp(logweights)));
end
grid = cell(1,L); [grid{:}] = ndgrid(1:O.NNodes);
O.QuadratureNodes = zeros(O.NNodes^L,L);
O.QuadratureLogWeights = zeros(O.NNodes^L,1);
for l = 1:L
    O.QuadratureNodes(:,l) = nodes(grid{l}(:));
    O.QuadratureLogWeights = O.QuadratureLogWeights+logweights(grid{l}(:));
end
if strcmp(O.Quadrature,'hermite')
    O.QuadratureLogWeights = O.QuadratureLogWeights-log(sum(exp(O.QuadratureLogWeights)));
end
O.QuadratureWeights = exp(O.QuadratureLogWeights);
O.QuadratureMass = sum(O.QuadratureWeights);
O.NLatent = L; O.NCategories = K; O.ItemFactor = factor; O.Anchor = anchor;
O.NObs = n; O.NP = n;
O.NVarA = size(I.Xa,2); O.NVarStr = size(I.Xstr,2); O.NVarLatent = size(I.Xlatent,2);
namesA = names(O,'NamesA',O.NVarA,'Xa');
namesStr = names(O,'NamesStr',O.NVarStr,'Xstr');
namesLat = names(O,'NamesLatent',O.NVarLatent,'Xlatent');
namesMea = names(O,'NamesMea',J,'item');
b = zeros(0,1); labels = cell(0,1); layout = struct();
layout.beta = (1:O.NVarA)'; b = zeros(O.NVarA,1);
for j = 1:O.NVarA, labels{end+1,1} = ['count:' namesA{j}]; end
cons = find(all(I.Xa==I.Xa(1,:),1) & I.Xa(1,:)~=0,1);
if ~isempty(cons), b(cons) = log(max(sum(I.W.*I.Y)/sum(I.W),0.1))/I.Xa(1,cons); end
layout.logalpha = [];
if O.NB, layout.logalpha = numel(b)+1; b(end+1,1) = log(0.5); labels{end+1,1} = 'log(alpha)'; end
layout.delta = zeros(0,L);
if O.LatentActive
    layout.delta = reshape(numel(b)+(1:O.NVarLatent*L),O.NVarLatent,L);
    b = [b;zeros(O.NVarLatent*L,1)];
    for l = 1:L
        for j = 1:O.NVarLatent, labels{end+1,1} = sprintf('LV%d*%s',l,namesLat{j}); end
    end
end
layout.gamma = reshape(numel(b)+(1:O.NVarStr*L),O.NVarStr,L);
b = [b;zeros(O.NVarStr*L,1)];
for l = 1:L
    for j = 1:O.NVarStr, labels{end+1,1} = sprintf('LV%d<- %s',l,namesStr{j}); end
end
layout.loading = (numel(b)+(1:J))';
b = [b;0.7*ones(J,1)]; b(layout.loading(anchor)) = log(0.7);
for j = 1:J
    prefix = 'loading'; if anchor(j), prefix = 'log(anchor loading)'; end
    labels{end+1,1} = sprintf('%s:%s',prefix,namesMea{j});
end
layout.thresholds = cell(J,1);
for j = 1:J
    layout.thresholds{j} = (numel(b)+(1:K(j)-1))';
    observed = ~isnan(I.Xmea(:,j)); cuts = zeros(K(j)-1,1);
    for k = 1:K(j)-1
        share = sum(I.W(observed).*(I.Xmea(observed,j)<=k))/sum(I.W(observed));
        cuts(k) = log(share/(1-share));
        prefix = 'threshold1'; if k>1, prefix = sprintf('log(gap%d)',k); end
        labels{end+1,1} = sprintf('%s:%s',prefix,namesMea{j});
    end
    b = [b;cuts(1);log(diff(cuts))];
end
O.Layout = layout; O.Params = numel(b); O.ParameterNames = labels;
O.HybridVersion = 1;
O.InitialValues = b;
O.NamesA = namesA; O.NamesStr = namesStr;
O.NamesLatent = namesLat; O.NamesMea = namesMea;
O.NVarMea = J; O.NVarNB = double(O.NB); O.NVarZinf = 0;
O.NCTMiss = ones(n,1); O.NumGrad = 0;
Results = struct();
OptimOpt = optimset('Display','off','GradObj','on','MaxIter',1000, ...
    'MaxFunEvals',10000,'TolFun',1e-7,'TolX',1e-7);
for field = {'SaveXLS','SaveTxtOutput'}
    if ~isfield(O,field{1}) || isempty(O.(field{1})), O.(field{1}) = double(O.Display~=0); end
    value = O.(field{1});
    if ~isscalar(value) || ~ismember(value,[0,1])
        error('CDM_hybrid:Options','%s must be 0 or 1.',field{1});
    end
end
if ~O.EvaluateOnly && (O.SaveXLS || O.SaveTxtOutput)
    if exist('setupDceOutputDefaults','file')~=2
        error('CDM_hybrid:MissingTools','Add the companion Tools folder for standard output, or set SaveXLS=0 and SaveTxtOutput=0.');
    end
    O = setupDceOutputDefaults(O);
end
end

function out = names(O,field,k,prefix)
if isfield(O,field) && ~isempty(O.(field))
    out = cellstr(string(O.(field)(:)));
    if numel(out)~=k, error('CDM_hybrid:Options','%s must contain %d labels.',field,k); end
else
    out = arrayfun(@(j) sprintf('%s%d',prefix,j),(1:k)','UniformOutput',false);
end
end

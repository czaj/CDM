function Results = CDM_hybrid_output(Results,exportOpt)
% CDM_hybrid_output Standard CDM/DCE screen and XLSX output for a hybrid fit.
%   R = CDM_hybrid_output(R,exportOpt) also exports an existing saved fit.
%   Output controls: OutputDir, OutputFilePrefix, ProjectName, SourceScript,
%   OutputLogFile, Display, xlsOverwrite and OutputNotes (an N-by-2 cell array).
%   Explicit calls write XLSX even with Display=0. No model is re-estimated.
%   Public SEs require a valid Hessian and an explicit successful quadrature
%   check. An existing InferenceValid=false or BoundaryCompetitive=true also
%   withholds inference. Local covariance diagnostics remain in the MAT fit.
%   Positive loadings, alpha and ordered cutoffs use the delta method.
%   See also CDM_hybrid, genOutput, writeResultsXlsx.

if nargin<2 || isempty(exportOpt), exportOpt = struct(); end
if ~isstruct(Results) || ~all(isfield(Results,{'EstimOpt','INPUT','bhat','Layout','ParameterNames','LL'}))
    error('CDM_hybrid_output:Input','Supply a CDM_hybrid result with saved inputs and parameter layout.');
end
O = Results.EstimOpt;
for name = {'OutputDir','OutputFilePrefix','ProjectName','SourceScript','OutputLogFile','Display','xlsOverwrite','OutputNotes'}
    key = name{1};
    if isfield(exportOpt,key), O.(key) = exportOpt.(key); end
end
if ~isfield(O,'Display'), O.Display = 1; end
if ~isfield(O,'OutputDir') || isempty(O.OutputDir)
    % setupDceOutputDefaults owns the estimator diary. Standalone export only
    % chooses its destination and must not turn off or replace that diary.
    if isfield(O,'SourceScript') && exist(char(O.SourceScript),'file')==2
        O.OutputDir = fullfile(fileparts(char(O.SourceScript)),'output');
    else
        O.OutputDir = fullfile(pwd,'output');
    end
end
if ~isfield(O,'xlsOverwrite'), O.xlsOverwrite = 1; end
if isfield(O,'OutputNotes') && (~iscell(O.OutputNotes) || size(O.OutputNotes,2)~=2)
    error('CDM_hybrid_output:Options','OutputNotes must be an N-by-2 cell array.');
end

b = Results.bhat(:); p = numel(b); L = Results.Layout;
V = nan(p);
if isfield(Results,'ihess') && isequal(size(Results.ihess),[p,p]), V = Results.ihess; end
hValid = isfield(Results,'HessDiagnostics') && isfield(Results.HessDiagnostics,'ValidInference') && ...
    isequal(Results.HessDiagnostics.ValidInference,true);
qChecked = flag(Results,'QuadratureChecked') || flag(Results,'QuadratureCheckPerformed');
qStable = flag(Results,'QuadratureStable');
valid = hValid && qChecked && qStable && all(isfinite(V(:)));
if isfield(Results,'InferenceValid'), valid = valid && flag(Results,'InferenceValid'); end
if flag(Results,'BoundaryCompetitive'), valid = false; end
Results.OutputInferenceValid = valid;
Results.Dist = -ones(numel(L.beta),1);
Names = struct(); Heads = struct(); template = cell(0,1);
O.NoSigRows = struct();
rawNames = cellstr(string(Results.ParameterNames(:)));
add('DetailsA',strrep(rawNames(L.beta),'count:',''),b(L.beta),eyeRows(L.beta), ...
    'Count equation',[]);
if ~isempty(L.logalpha)
    a = exp(b(L.logalpha)); D = zeros(1,p); D(L.logalpha) = a;
    add('DetailsNB',{'alpha'},a,D,'NB2 overdispersion (natural scale)',Inf);
end
for l = 1:size(L.delta,2)
    idx = L.delta(:,l);
    if ~isempty(idx)
        add(sprintf('DetailsLV%d',l),rawNames(idx),b(idx),eyeRows(idx), ...
            sprintf('LV%d effects in the count equation',l),[]);
    end
end
for l = 1:size(L.gamma,2)
    idx = L.gamma(:,l);
    if ~isempty(idx)
        add(sprintf('DetailsSE%d',l),rawNames(idx),b(idx),eyeRows(idx), ...
            sprintf('LV%d structural equation',l),[]);
    end
end
itemNames = labels(O,'NamesMea',numel(L.loading),'item');
for j = 1:numel(L.loading)
    idx = L.thresholds{j}; nk = numel(idx);
    theta = [b(L.loading(j)); b(idx(1)); b(idx(1))+cumsum(exp(b(idx(2:end))))];
    D = zeros(nk+1,p); D(1,L.loading(j)) = 1;
    if O.Anchor(j)
        theta(1) = exp(theta(1)); D(1,L.loading(j)) = theta(1);
    end
    D(2:end,idx(1)) = 1;
    for k = 2:nk, D(k+1,idx(2:k)) = exp(b(idx(2:k))); end
    rowNames = [{sprintf('LV%d loading',O.ItemFactor(j))}; ...
        arrayfun(@(k) sprintf('Cutoff %d',k),(1:nk)','UniformOutput',false)];
    suppress = []; if O.Anchor(j), suppress = 1; end
    add(sprintf('DetailsMea%d',j),rowNames,theta,D, ...
        ['Measurement: ' itemNames{j}],suppress);
end
% The raw block makes optimizer starts, all log gaps and transformations
% auditable; its zero tests are omitted because several are not substantive.
add('DetailsRaw',rawNames,b,eye(p),'Raw optimizer parameters (no zero tests)',Inf);
if ~valid && any(isfinite(diag(V)))
    Results.DetailsLocalUnchecked = [b,zeros(p,1),sqrt(max(diag(V),0)),nan(p,1)];
    Names.DetailsLocalUnchecked = rawNames;
    Heads.DetailsLocalUnchecked = {'Local SE only - inference not validated';'tb'};
    template{end+1,1} = 'DetailsLocalUnchecked';
    O.NoSigRows.DetailsLocalUnchecked = Inf;
end

n = numel(Results.INPUT.Y);
AIC = 2*p-2*Results.LL; BIC = log(n)*p-2*Results.LL;
LL0 = NaN; if isfield(Results,'LL0') && isscalar(Results.LL0), LL0 = Results.LL0; end
Results.stats = [Results.LL;LL0;1-Results.LL/LL0;NaN;AIC/n;BIC/n;n;n;p];
Tail = cell(16,2);
Tail(2:11,1) = {'Model diagnostics';'LL at convergence';'LL at constant(s) only'; ...
    ['McFadden''s pseudo-R' char(178)];['Ben-Akiva-Lerman''s pseudo-R' char(178)]; ...
    'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)'};
Tail(3:11,2) = num2cell(Results.stats);
Tail(13:16,:) = {'Estimation method','Joint maximum likelihood'; ...
    'Optimization method','fminunc';'Gradient','Analytic';'Hessian','Numerical derivative of analytic score'};
latentLaw = 'population'; if isfield(O,'LatentDistribution'), latentLaw = O.LatentDistribution; end
countIndex = 'log(mu)'; if flag(O,'LogSeries'), countIndex = 'log(alpha * mu), NB2 boundary'; end
covType = 'Observed information'; if flag(O,'RobustStd'), covType = 'Observation-score sandwich'; end
Tail = [Tail; {'Joint AIC',AIC;'Joint BIC',BIC;'Count index',countIndex; ...
    'Normal factor law',latentLaw;'Positive-count sample',yesno(flag(O,'Truncated')); ...
    'Latent count link',yesno(flag(O,'LatentActive'));'Covariance',covType; ...
    'Hessian / convergence valid',yesno(hValid);'Quadrature check performed',yesno(qChecked); ...
    'Quadrature stable',yesno(qStable);'Public inference validated',yesno(valid); ...
    'Finite-NB2 boundary competitive',yesno(flag(Results,'BoundaryCompetitive')); ...
    'Integration rule',O.Quadrature;'Nodes per factor',num2str(O.NNodes); ...
    'Integration points',num2str(numel(O.QuadratureWeights)); ...
    'AIC/BIC comparisons','Same rows and same indicator model'; ...
    'Cutoff p-values','Diagnostic zero tests only'; ...
    'alpha / positive anchor','No zero tests or significance stars'; ...
    'Structural intercepts','Fixed at zero for identification'; ...
    'Latent residual covariance','Identity (factor scale identification)'}];
if strcmp(O.Quadrature,'normal-grid'), Tail(end+1,:) = {'Grid range',sprintf('[-%g, %g]',O.GridLimit,O.GridLimit)}; end
if isfield(Results,'exitflag'), Tail(end+1,:) = {'Optimizer exit flag',num2str(Results.exitflag)}; end
if isfield(Results,'HessDiagnostics')
    H = Results.HessDiagnostics;
    if isfield(H,'ReciprocalCondition'), Tail(end+1,:) = {'Hessian reciprocal condition',sprintf('%.5g',H.ReciprocalCondition)}; end
    if isfield(H,'GradientPerWeight'), Tail(end+1,:) = {'Max gradient / weight',sprintf('%.5g',H.GradientPerWeight)}; end
end
if isfield(O,'OutputNotes'), Tail = [Tail;cell(1,2);wrapNotes(O.OutputNotes)]; end
% All footer values after the diagnostics are printed as strings by genOutput.
for k = 13:size(Tail,1)
    for c = 1:2
        if isempty(Tail{k,c}), Tail{k,c} = ''; elseif isnumeric(Tail{k,c}), Tail{k,c} = num2str(Tail{k,c},12); end
    end
end
subtitle = sprintf('%s-normal factors; %d latent factor(s)',latentLaw,O.NLatent);
Head = {Results.Model,subtitle};
Results.EstimOpt = O;
if isfield(O,'OutputLogFile'), Results.output_txt = O.OutputLogFile; end
if ~isfield(Results,'EstimationSeconds'), Results.EstimationSeconds = NaN; end
[Results.R_out,Results.OutputXLSX] = genOutput(O,Results,Head,Tail,Names,template,template,Heads,template);
Results.OutputNames = Names;
Results.OutputBlocks = template;

    function D = eyeRows(idx)
        D = zeros(numel(idx),p);
        D(sub2ind(size(D),(1:numel(idx))',idx(:))) = 1;
    end

    function add(key,rowNames,value,D,heading,noSig)
        se = nan(numel(value),1);
        if valid, se = sqrt(max(sum((D*V).*D,2),0)); end
        pv = erfc(abs(value(:)./se)/sqrt(2));
        Results.(key) = [value(:),zeros(numel(value),1),se,pv];
        if ~isempty(noSig)
            O.NoSigRows.(key) = noSig;
            if any(isinf(noSig)), Results.(key)(:,4) = NaN; else, Results.(key)(noSig,4) = NaN; end
        end
        Names.(key) = rowNames(:); Heads.(key) = {heading;'tb'};
        template{end+1,1} = key;
    end
end

function tf = flag(S,key)
tf = isfield(S,key) && isequal(S.(key),true);
end

function out = yesno(tf)
out = 'No'; if tf, out = 'Yes'; end
end

function out = labels(O,key,n,prefix)
if isfield(O,key) && ~isempty(O.(key)), out = cellstr(string(O.(key)(:)));
else, out = arrayfun(@(j) sprintf('%s%d',prefix,j),(1:n)','UniformOutput',false); end
end

function out = wrapNotes(notes)
% Keep the shared Excel AutoFit from turning a source path into a huge column.
out = cell(0,2);
for i = 1:size(notes,1)
    label = notes{i,1}; value = char(string(notes{i,2}));
    while numel(value)>64
        cut = find(isspace(value(1:64)),1,'last');
        if isempty(cut) || cut<32, cut=64; end
        out(end+1,:) = {label,strtrim(value(1:cut))}; %#ok<AGROW>
        label=''; value=strtrim(value(cut+1:end));
    end
    out(end+1,:) = {label,value}; %#ok<AGROW>
end
end

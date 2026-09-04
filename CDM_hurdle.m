function Results = CDM_hurdle(INPUT,Results_old,EstimOpt,OptimOpt)
% CDM_hurdle estimates a two-part (hurdle) count data model.
%
% Syntax:   Results = CDM_hurdle(INPUT,Results_old,EstimOpt,OptimOpt)
%
% The first part is a binary model of whether the respondent participates at all;
% the second is a zero-truncated count model for how many trips a participant
% takes. For recreation demand this is usually better motivated than zero
% inflation: not going at all is a different decision from deciding how often to
% go, rather than a hidden class.
%
% The log-likelihood separates,
%
%   LL = sum_{y=0} log(1-P) + sum_{y>0} [ log(P) + log f(y | y>0) ]
%
% and the two parts share no parameters, so each is maximised on its own and the
% information matrix is block diagonal. That is not a shortcut - it is what the
% model implies.
%
% INPUT:
%    Y, Xa, Xnb, W, MissingInd - as for CDM
%    Xh - covariates of participation; a constant is prepended. Defaults to Xa,
%         which is the usual specification
%
% EstimOpt:
%    NB, NBtype, Censored, weights and the estimation options behave as in CDM
%    and apply to the count part. Truncated is forced to 1 - that is the model.
%    probit - 1 for a probit participation equation, otherwise logit
%    NamesH - names of the participation covariates (without the constant)
%
% The result is shaped like a CDM result, with the participation equation stored
% in the zero-inflation slot (P = 1-q), so CDM_post and CDM_compare work on it
% unchanged. Expected trips are P*E[Y|Y>0], which is what CDM_EY already returns
% for this combination.
%
% See also CDM, CDM_post, CDM_compare, CDM_EY
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

if nargin < 4
    error('Too few input arguments for CDM_hurdle(INPUT,Results_old,EstimOpt,OptimOpt)')
end
if ~isfield(EstimOpt,'Display')
    EstimOpt.Display = 1;
end

Y = INPUT.Y(:);
d = double(Y > 0); % participation indicator
n = numel(Y);

% Same design screen as CDM - the count part of a hurdle model reads its
% consumer surplus from one coefficient in exactly the same way.
DesignCheck = CDM_checkdesign(INPUT,EstimOpt);

if isfield(INPUT,'Xh') && ~isempty(INPUT.Xh)
    Xh = INPUT.Xh;
else
    Xh = INPUT.Xa;
end
isConst = all(Xh == Xh(1,:),1);
Xh(:,isConst) = [];
Xh = [ones(n,1),Xh]; % the package owns the constant, as it does for Xnb and Xzinf
NVarH = size(Xh,2);

if isfield(INPUT,'W') && ~isempty(INPUT.W) && numel(INPUT.W) == n
    W = INPUT.W(:);
else
    W = ones(n,1);
end

useProbit = isfield(EstimOpt,'probit') && EstimOpt.probit == 1;

if EstimOpt.Display ~= 0
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
    if useProbit
        disp('Estimating hurdle model: probit participation ...')
    else
        disp('Estimating hurdle model: logit participation ...')
    end
end

%% Part one - participation

Opt = optimoptions('fminunc','Algorithm','quasi-newton','GradObj','on', ...
    'Display','off','MaxIter',1e4,'MaxFunEvals',1e5);
dh0 = zeros(NVarH,1);
dh0(1) = log(max(mean(d),1e-3)/max(1-mean(d),1e-3)); % start at the observed share

[dhat,nLLh,exitH,~,~,hessH] = fminunc(@(D) binaryLL(D,d,Xh,W,useProbit),dh0,Opt);

Results.DesignCheck = DesignCheck;
Results.Binary = struct();
Results.Binary.bhat = dhat;
Results.Binary.LL = -nLLh;
Results.Binary.exitf = exitH;
[Results.Binary.ihess,Results.Binary.HessDiagnostics] = hessianInverse(hessH,'hurdle participation');
Results.Binary.std = sqrt(diag(Results.Binary.ihess));
Results.Binary.probit = useProbit;

%% Part two - trips, conditional on participating

INPUT2 = INPUT;
keep = d == 1;
flds = fieldnames(INPUT2);
for i = 1:numel(flds)
    tmp = INPUT2.(flds{i});
    if size(tmp,1) == n
        INPUT2.(flds{i}) = tmp(keep,:);
    end
end
if isfield(INPUT2,'Xh')
    INPUT2 = rmfield(INPUT2,'Xh');
end

EstimOpt2 = EstimOpt;
EstimOpt2.Truncated = 1;   % this is the model
EstimOpt2.Zinf = 0;
EstimOpt2.NP = sum(keep);
EstimOpt2.NCT = 1;
EstimOpt2.NObs = sum(keep);
EstimOpt2.NCTMiss = ones(sum(keep),1);
if isfield(INPUT2,'MissingInd')
    INPUT2.MissingInd = zeros(sum(keep),1);
end

Results.Count = CDM(INPUT2,Results_old,EstimOpt2,OptimOpt);

%% Assemble, in the shape of a CDM result

NVarA = Results.Count.EstimOpt.NVarA;
NVarNB = Results.Count.EstimOpt.NVarNB;

% P = 1-q, so storing -dhat in the zero-inflation slot makes q = 1-P
Results.bhat = [Results.Count.bhat; -dhat];
kCount = numel(Results.Count.bhat);
Results.ihess = zeros(kCount+NVarH);
Results.ihess(1:kCount,1:kCount) = Results.Count.ihess;
Results.ihess(kCount+1:end,kCount+1:end) = Results.Binary.ihess; % block diagonal by construction
Results.std = sqrt(diag(Results.ihess));
Results.R = [Results.bhat,Results.std,pv(Results.bhat,Results.std)];

Results.LL = Results.Binary.LL + Results.Count.LL;
Results.LLdetailed = zeros(n,1);
Ph = participationProb(dhat,Xh,useProbit);
Results.LLdetailed(~keep) = W(~keep).*log(max(1-Ph(~keep),realmin));
Results.LLdetailed(keep) = W(keep).*log(max(Ph(keep),realmin)) + Results.Count.LLdetailed;

EstimOpt.Truncated = 1;
EstimOpt.Zinf = 1;      % the participation equation lives in the inflation slot
EstimOpt.Hurdle = 1;
EstimOpt.NB = Results.Count.EstimOpt.NB;
EstimOpt.NBtype = Results.Count.EstimOpt.NBtype;
EstimOpt.Censored = Results.Count.EstimOpt.Censored;
EstimOpt.RealMin = Results.Count.EstimOpt.RealMin;
EstimOpt.NVarA = NVarA;
EstimOpt.NVarNB = NVarNB;
EstimOpt.NVarZinf = NVarH;
EstimOpt.NObs = n;
EstimOpt.Params = kCount + NVarH;
EstimOpt.probit = useProbit;
EstimOpt.NamesA = Results.Count.EstimOpt.NamesA;
EstimOpt.NamesNB = Results.Count.EstimOpt.NamesNB;
EstimOpt.TC_par_index = Results.Count.EstimOpt.TC_par_index;

if isfield(EstimOpt,'NamesH') && ~isempty(EstimOpt.NamesH) && numel(EstimOpt.NamesH) == NVarH-1
    EstimOpt.NamesZinf = [{'Cons'};EstimOpt.NamesH(:)];
else
    EstimOpt.NamesZinf = [{'Cons'};cellstr(num2str((1:NVarH-1)'))];
end

INPUT.Xzinf = Xh;
INPUT.Xnb = Results.Count.INPUT.Xnb;
if size(INPUT.Xnb,1) ~= n % the count part only saw participants
    if NVarNB > 0
        INPUT.Xnb = ones(n,1);
    else
        INPUT.Xnb = zeros(n,0);
    end
end
INPUT.W = W;

Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;

% constants-only hurdle, for McFadden's pseudo-R2
share = sum(W.*d)/sum(W);
LL0h = sum(W.*(d*log(max(share,realmin)) + (1-d)*log(max(1-share,realmin))));
if isfield(Results_old,'CDM0') && isfield(Results_old.CDM0,'LL')
    LL0 = LL0h + Results_old.CDM0.LL;
    k = EstimOpt.Params;
    Results.stats = [Results.LL; LL0; 1-Results.LL/LL0; NaN; ...
        (2*k-2*Results.LL)/n; (log(n)*k-2*Results.LL)/n; n; EstimOpt.NP; k];
else
    Results.stats = [];
end

%% Print

if EstimOpt.Display ~= 0
    fprintf('\n')
    cprintf('*Black','Hurdle model \n')
    fprintf('%-30s%14s%12s%10s\n','','coef.','st.err.','p-value')
    cprintf('*Black','Participation \n')
    for i = 1:NVarH
        fprintf('%-30s%14.4f%12.4f%10.4f\n',char(EstimOpt.NamesZinf(i)), ...
            dhat(i),Results.Binary.std(i),pv(dhat(i),Results.Binary.std(i)))
    end
    cprintf('*Black','Trips, given participation \n')
    for i = 1:NVarA
        fprintf('%-30s%14.4f%12.4f%10.4f\n',char(EstimOpt.NamesA(i)), ...
            Results.Count.bhat(i),Results.Count.std(i),pv(Results.Count.bhat(i),Results.Count.std(i)))
    end
    for i = 1:NVarNB
        fprintf('%-30s%14.4f%12.4f%10s\n',['alpha: ' char(EstimOpt.NamesNB(i))], ...
            Results.Count.bhat(NVarA+i),Results.Count.std(NVarA+i),'-')
    end
    fprintf('\n')
    fprintf('%-30s%14.4f\n','LL participation:',Results.Binary.LL)
    fprintf('%-30s%14.4f\n','LL trips:',Results.Count.LL)
    fprintf('%-30s%14.4f\n','LL total:',Results.LL)
    if ~isempty(Results.stats)
        fprintf('%-30s%14.4f\n','McFadden pseudo-R2:',Results.stats(3))
        fprintf('%-30s%14.4f\n','AIC/n:',Results.stats(5))
    end
    fprintf('%-30s%14d\n','k (parameters):',EstimOpt.Params)
    fprintf('\n')
end
end


function P = participationProb(D,Xh,useProbit)
if useProbit
    P = normcdf(Xh*D);
else
    P = 1./(1+exp(-Xh*D));
end
P = min(max(P,realmin),1-eps);
end


function [nLL,g] = binaryLL(D,d,Xh,W,useProbit)
% Negative weighted log-likelihood of the participation equation, with gradient.
P = participationProb(D,Xh,useProbit);
nLL = -sum(W.*(d.*log(P) + (1-d).*log(1-P)));
if nargout > 1
    if useProbit
        dens = normpdf(Xh*D);
        res = dens.*(d-P)./(P.*(1-P));
    else
        res = d-P;
    end
    g = -(Xh'*(W.*res));
end
end

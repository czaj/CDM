function Results = CDM_post(Results,EstimOpt)
% CDM_post post-estimation results for a count data model: consumer surplus,
% average marginal effects, elasticities and the fit to the count distribution.
%
% Syntax:   Results = CDM_post(Results)
%           Results = CDM_post(Results,EstimOpt)
%
% Results must come from CDM. Options are taken from Results.EstimOpt and can be
% overridden by the second argument:
%
%   NSim        number of Krinsky-Robb draws (default 1e3)
%   CIlevel     confidence level for the simulated intervals (default 0.95)
%   Display     0 to compute quietly
%   TC_par_index  index of the travel cost parameter; consumer surplus is skipped
%                 when it is not set
%
% Adds to Results:
%   EY          expected trips per observation, at the estimates
%   CS          consumer surplus per trip and per person, with delta-method and
%               Krinsky-Robb intervals, plus the sample aggregate
%   Margins     average marginal effects on E[Y] and elasticities, with
%               Krinsky-Robb standard errors
%   CountFit    observed against predicted count frequencies, Pearson chi-square
%               and the deviance pseudo-R2
%
% Consumer surplus per trip is -1/beta_tc (Hellerstein and Mendelsohn 1993). Per
% person it is E[Y]/(-beta_tc), where E[Y] follows the sampling scheme - see
% CDM_EY. Krinsky-Robb intervals are reported alongside the delta method because
% the delta method is unreliable for a ratio whose denominator is estimated.
%
% See also CDM, CDM_EY, CDM_compare
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

if nargin < 2 || isempty(EstimOpt)
    EstimOpt = struct();
end
EO = Results.EstimOpt;
fn = fieldnames(EstimOpt);
for i = 1:numel(fn)
    EO.(fn{i}) = EstimOpt.(fn{i});
end
if ~isfield(EO,'NSim') || isempty(EO.NSim)
    EO.NSim = 1e3;
end
if ~isfield(EO,'CIlevel') || isempty(EO.CIlevel)
    EO.CIlevel = 0.95;
end
if ~isfield(EO,'Display')
    EO.Display = 1;
end
alpha = (1-EO.CIlevel)/2;

IN = Results.INPUT;
b = Results.bhat(:);
W = IN.W(:);
W = W/sum(W); % weights that average to one observation

if ~isfield(IN,'Xnb'),   IN.Xnb   = zeros(size(IN.Y,1),0); end
if ~isfield(IN,'Xzinf'), IN.Xzinf = zeros(size(IN.Y,1),0); end

%% Krinsky-Robb draws of the parameter vector

rng(EO.NSim + 20260829,'twister'); % reproducible across runs
V = Results.ihess;
V(~isfinite(V)) = 0;
V = (V+V')/2;
[L,p] = chol(V,'lower');
if p ~= 0
    L = chol(nearestSPD(V),'lower');
end
Bsim = b + L*randn(numel(b),EO.NSim);

%% Expected trips and consumer surplus

Results.EY = CDM_EY(b,IN.Xa,IN.Xnb,IN.Xzinf,EO);

TCidx = 0;
if isfield(EO,'TC_par_index') && ~isempty(EO.TC_par_index)
    TCidx = EO.TC_par_index;
end

EYsim = zeros(EO.NSim,1);
CSsim = NaN(EO.NSim,2);
for r = 1:EO.NSim
    EYr = CDM_EY(Bsim(:,r),IN.Xa,IN.Xnb,IN.Xzinf,EO);
    EYsim(r) = W'*EYr;
    if TCidx > 0
        csTrip = -1/Bsim(TCidx,r);
        CSsim(r,:) = [csTrip, EYsim(r)*csTrip];
    end
end

Results.CS = struct();
Results.CS.EYmean = W'*Results.EY;
Results.CS.EYsim = quantile(EYsim,[alpha 0.5 1-alpha]);

if TCidx > 0 && TCidx <= EO.NVarA
    bTC = b(TCidx);
    seTC = sqrt(Results.ihess(TCidx,TCidx));
    csTrip = -1/bTC;
    csTripSE = sqrt(Results.ihess(TCidx,TCidx)/bTC^4); % delta method
    csPers = Results.CS.EYmean*csTrip;

    Results.CS.PerTrip = [csTrip, csTripSE, quantile(CSsim(:,1),[alpha 0.5 1-alpha])];
    Results.CS.PerPerson = [csPers, NaN, quantile(CSsim(:,2),[alpha 0.5 1-alpha])];
    % Inverting the confidence interval on beta is the exact interval for a ratio
    % with a constant numerator, i.e. what Fieller reduces to in this case.
    bCI = bTC + norminv([1-alpha alpha])*seTC;
    Results.CS.PerTripFieller = sort(-1./bCI); % -1/beta rises with beta, so a negative beta flips the ends
    Results.CS.Aggregate = sum(IN.W(:))*csPers/numel(IN.W); % per respondent, weighted
    Results.CS.TCname = char(EO.NamesA(TCidx));
end

%% Average marginal effects and elasticities on E[Y]

isConst = all(IN.Xa == IN.Xa(1,:),1);
okElast = all(IN.Xa > 0,1); % an elasticity is only meaningful for a strictly positive regressor
NVarA = EO.NVarA;
AMEsim = zeros(EO.NSim,NVarA);
ELAsim = zeros(EO.NSim,NVarA);
step = 1e-5*max(1,std(IN.Xa,0,1));

for r = 1:EO.NSim
    br = Bsim(:,r);
    EYr = CDM_EY(br,IN.Xa,IN.Xnb,IN.Xzinf,EO);
    for j = 1:NVarA
        if isConst(j)
            AMEsim(r,j) = NaN; ELAsim(r,j) = NaN;
            continue
        end
        Xp = IN.Xa; Xp(:,j) = Xp(:,j)+step(j);
        Xm = IN.Xa; Xm(:,j) = Xm(:,j)-step(j);
        d = (CDM_EY(br,Xp,IN.Xnb,IN.Xzinf,EO)-CDM_EY(br,Xm,IN.Xnb,IN.Xzinf,EO))/(2*step(j));
        AMEsim(r,j) = W'*d;
        if okElast(j)
            ELAsim(r,j) = W'*(d.*IN.Xa(:,j)./EYr);
        else
            ELAsim(r,j) = NaN;
        end
    end
end

AME = zeros(NVarA,1); ELA = zeros(NVarA,1);
EYb = Results.EY;
for j = 1:NVarA
    if isConst(j)
        AME(j) = NaN; ELA(j) = NaN;
        continue
    end
    Xp = IN.Xa; Xp(:,j) = Xp(:,j)+step(j);
    Xm = IN.Xa; Xm(:,j) = Xm(:,j)-step(j);
    d = (CDM_EY(b,Xp,IN.Xnb,IN.Xzinf,EO)-CDM_EY(b,Xm,IN.Xnb,IN.Xzinf,EO))/(2*step(j));
    AME(j) = W'*d;
    if okElast(j)
        ELA(j) = W'*(d.*IN.Xa(:,j)./EYb);
    else
        ELA(j) = NaN;
    end
end

Results.Margins = struct('AME',[AME,std(AMEsim,0,1)'], ...
                         'Elasticity',[ELA,std(ELAsim,0,1)'], ...
                         'Names',{EO.NamesA});

%% Fit to the count distribution

% P(Y = k) comes straight from the estimated likelihood, so every variant -
% truncation, stratification, censoring, zero inflation - is handled consistently.
kMin = double(EO.Truncated > 0);
kMax = max(IN.Y);
kk = kMin:kMax;
Pk = zeros(numel(kk),1);
for i = 1:numel(kk)
    Pk(i) = W'*exp(-LL_CDM(kk(i)*ones(size(IN.Y)),IN.Xa,IN.Xnb,IN.Xzinf,EO,b));
end
nObs = numel(IN.Y);
Obs = zeros(numel(kk),1);
for i = 1:numel(kk)
    Obs(i) = sum(IN.W(IN.Y == kk(i)))/sum(IN.W);
end
Exp = Pk;
keep = Exp*nObs >= 1; % cells too thin for a chi-square contribution
chi2 = nObs*sum((Obs(keep)-Exp(keep)).^2./Exp(keep));
df = max(1,sum(keep)-1-EO.Params);

Results.CountFit = struct('k',kk(:),'Observed',Obs,'Predicted',Exp, ...
    'Chi2',chi2,'df',df,'pValue',1-chi2cdf(chi2,df));

%% Print

if EO.Display ~= 0
    fprintf('\n')
    cprintf('*Black','Post-estimation results \n')
    fprintf('%-30s%12.4f\n','Expected trips E[Y]:',Results.CS.EYmean)
    fprintf('%-30s%22s\n','  simulated 95% interval:', ...
        sprintf('[%.4f, %.4f]',Results.CS.EYsim(1),Results.CS.EYsim(3)))

    if isfield(Results.CS,'PerTrip')
        fprintf('\n')
        cprintf('*Black','Consumer surplus \n')
        fprintf('%-30s%12s\n','travel cost variable:',Results.CS.TCname)
        fprintf('%-16s%12s%12s%22s%22s\n','','value','s.e.','Krinsky-Robb 95%','inverted CI on beta')
        fprintf('%-16s%12.4f%12.4f%22s%22s\n','per trip', ...
            Results.CS.PerTrip(1),Results.CS.PerTrip(2), ...
            sprintf('[%.3f, %.3f]',Results.CS.PerTrip(3),Results.CS.PerTrip(5)), ...
            sprintf('[%.3f, %.3f]',Results.CS.PerTripFieller(1),Results.CS.PerTripFieller(2)))
        fprintf('%-16s%12.4f%12s%22s\n','per person', ...
            Results.CS.PerPerson(1),'-', ...
            sprintf('[%.3f, %.3f]',Results.CS.PerPerson(3),Results.CS.PerPerson(5)))
    else
        fprintf('\n')
        cprintf(rgb('DarkOrange'),'Consumer surplus not reported - set EstimOpt.TC_par_index \n')
    end

    fprintf('\n')
    cprintf('*Black','Average marginal effects on E[Y] \n')
    fprintf('%-30s%12s%10s%12s%10s\n','var.','dE[Y]/dx','s.e.','elasticity','s.e.')
    for j = 1:NVarA
        if isnan(AME(j))
            continue
        end
        if isnan(ELA(j))
            fprintf('%-30s%12.4f%10.4f%12s%10s\n',char(EO.NamesA(j)),AME(j),Results.Margins.AME(j,2),'-','-')
        else
            fprintf('%-30s%12.4f%10.4f%12.4f%10.4f\n',char(EO.NamesA(j)), ...
            AME(j),Results.Margins.AME(j,2),ELA(j),Results.Margins.Elasticity(j,2))
        end
    end

    fprintf('\n')
    cprintf('*Black','Fit to the count distribution \n')
    fprintf('%-8s%14s%14s\n','y','observed','predicted')
    for i = 1:min(numel(kk),11)
        fprintf('%-8d%14.4f%14.4f\n',kk(i),Obs(i),Exp(i))
    end
    if numel(kk) > 11
        fprintf('%-8s%14.4f%14.4f\n',['>' num2str(kk(11))],sum(Obs(12:end)),sum(Exp(12:end)))
    end
    fprintf('%-8s%14s%14.4f  (df %d, p %.4f)\n','Pearson','chi2',chi2,df,Results.CountFit.pValue)
    fprintf('\n')
end
end

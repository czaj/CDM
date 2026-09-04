function Results = CDM(INPUT,Results_old,EstimOpt,OptimOpt)
% CDM estimates a count data model - Poisson or Negative Binomial, optionally
% zero-inflated, zero-truncated (with or without endogenous stratification)
% and right-censored. Intended for single-site travel cost method studies.
%
% Syntax:   CDM(INPUT,EstimOpt,OptimOpt)
%           CDM(INPUT,Results_old,EstimOpt,OptimOpt)
%
% Inputs:
%    INPUT - clean, updated INPUT data from DataCleanCDM
%    EstimOpt - Estimation Options (see below)
%    OptimOpt - Optimizer Options define how the algorithm converges to the final
%    result. They are set by default based on the EstimOpt provided to
%    DataCleanCDM, however, they are subject to change.
%    Results_old - old results, used for starting values and for the LL of the
%    constants-only model (Results.CDM0, produced by DataCleanCDM)
%
% INPUT data:
%    Y - counts (non-negative integers), NCT*NP rows
%    Xa - explanatory variables of the count equation (include your own constant)
%    Xnb - optional covariates of overdispersion (a constant is prepended)
%    Xzinf - optional covariates of zero inflation (a constant is prepended)
%    W - optional weights, one per respondent or one per observation
%    MissingInd - optional 0/1 indicator; flagged observations are removed
%
% EstimOpt Options:
% Set them by e.g. EstimOpt.ProjectName = 'Project'
%
% General basics:
%   Display - 1; shows output, set to 0 to hide it
%   ProjectName - name of the project/model, appended to the saved .xls file
%   NCT - number of observations per respondent (1 for a cross-section)
%   NP - number of respondents
%
% Model options:
%   NB - 0 for Poisson (default), 1 for Negative Binomial
%   NBtype - 2 for NB2, variance lambda + alpha*lambda^2 (default), 1 for NB1,
%     variance lambda*(1+alpha), 3 for NB-P, variance lambda + alpha*lambda^P with
%     the power P estimated. NB1 and NB2 are NB-P at P = 1 and P = 2
%   Zinf - set to 1 for a zero-inflated model (requires INPUT.Xzinf)
%   probit - set to 1 to use a probit rather than a logit inflation equation
%   Truncated - 0 none (default), 1 zero truncation, 2 zero truncation with
%     endogenous stratification (on-site sample)
%   Censored - right-censoring point; observations above it are censored
%   TC_par_index - which parameter of Xa is the travel cost, used to report
%     consumer surplus per trip as -1/beta
%   AllowRankDeficient - set to 1 to downgrade the rank check on the explanatory
%     variables from an error to a warning (see CDM_checkdesign)
%
% Variables options:
%   NamesA - names of the variables in Xa, e.g. {'Cons';'Travel cost'}
%   NamesNB - names of the covariates of overdispersion (without the constant)
%   NamesZinf - names of the covariates of zero inflation (without the constant)
%
% Parameters options:
%   BActive = vector of 1s; set an element to 0 to constrain that parameter to
%     its initial value
%   ConstVarActive = 0; set to 1 to use constrained optimization
%
% Modelling options from DataCleanCDM:
%   ApproxHess = 1; for user supplied hessians, 1 for BHHH, 0 for analytical
%   RobustStd = 0; by default not using robust standard errors, set to 1 to use them
%   NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical gradient
%   RealMin = 0; set to 1 to guard the likelihood against underflow to zero
%   HessEstFix = 0; Options:
% o   0 - use optimization Hessian,
% o   1 - use jacobian-based (BHHH) Hessian,
% o   2 - use high-precision jacobian-based (BHHH) Hessian,
% o   3 - use numerical Hessian,
% o   4 - use analytical Hessian (not derived for count models - falls back to 3)
%
% Example:
%    [INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
%    Results.NB = CDM(INPUT,Results,EstimOpt,OptimOpt);
%
% See also DataCleanCDM, LL_CDM, CDM_demo
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

global B_backup;

%% 
tic 

Results.bhat = [];
Results.R = [];
%Results.R_out = {};
Results.stats = [];

%% check input:

if nargin<3
    error('Too few input arguments for CDM(INPUT,EstimOpt,OptimOpt)')
end
%how the output is presented
format shortG;
format compact;

if isfield(EstimOpt, 'NB') == 0 || EstimOpt.NB ==0
    EstimOpt.NB = 0;
    EstimOpt.NVarNB = 0;
end

if isfield(EstimOpt, 'Zinf') == 0 || EstimOpt.Zinf ==0
    EstimOpt.Zinf = 0;
    EstimOpt.NVarZinf = 0;
end

if isfield(EstimOpt,'TC_par_index') == 0 || isempty(EstimOpt.TC_par_index) || EstimOpt.TC_par_index == 0
    % Guessing index 2 used to be the default, which silently reported the consumer
    % surplus of whatever happened to sit there. Match on the name instead, and
    % report no consumer surplus at all rather than a confident wrong number.
    EstimOpt.TC_par_index = 0;
    if isfield(EstimOpt,'NamesA') && ~isempty(EstimOpt.NamesA)
        hit = find(~cellfun(@isempty,regexpi(EstimOpt.NamesA(:),'cost|koszt|price|fee')),1);
        if ~isempty(hit)
            EstimOpt.TC_par_index = hit;
        end
    end
end

TC_par_index = EstimOpt.TC_par_index;


% if EstimOpt.Zinf==1
%     INPUT.Xzinf = INPUT.Xzinf;
% else
%     INPUT.Xzinf = zeros(size(INPUT.Y,1),0);
% end

EstimOpt.NVarA = size(INPUT.Xa,2);

if isfield(INPUT, 'Xnb') == 0 || numel(INPUT.Xnb) == 0
    if EstimOpt.NB > 0
        INPUT.Xnb = ones(size(INPUT.Y,1),1);
    else
        INPUT.Xnb = zeros(size(INPUT.Y,1),0);
    end
else
    INPUT.Xnb = addConstant(INPUT.Xnb,'Xnb',EstimOpt.Display);
end
EstimOpt.NVarNB = size(INPUT.Xnb,2); % Covariates of alpha (overdispersion)

if EstimOpt.NVarNB > 0 && EstimOpt.NB == 0
   EstimOpt.NB = 1;
   disp('NVarNB are non zero. EstimOpt.NB = 1 has been set. ');
end

if isfield(INPUT, 'Xzinf') == 0 || numel(INPUT.Xzinf) == 0
    if EstimOpt.Zinf > 0
        INPUT.Xzinf = ones(size(INPUT.Y,1),1);
    else
        INPUT.Xzinf = zeros(size(INPUT.Y,1),0);
    end
else
    % Same convention as Xnb: supply the covariates, the constant is added here.
    INPUT.Xzinf = addConstant(INPUT.Xzinf,'Xzinf',EstimOpt.Display);
end

EstimOpt.NVarZinf = size(INPUT.Xzinf,2);

if EstimOpt.Display == 1
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
    if EstimOpt.NB == 0
        disp('Estimating Poisson regression model ...')
    elseif EstimOpt.NB == 1
        if isfield(EstimOpt,'NBtype') && EstimOpt.NBtype == 1
            disp('Estimating Negative Binomial (NB1, variance proportional to the mean) regression model ...')
        elseif isfield(EstimOpt,'NBtype') && EstimOpt.NBtype == 3
            disp('Estimating Negative Binomial (NB-P, power estimated) regression model ...')
        else
            disp('Estimating Negative Binomial (NB2) regression model ...')
        end
    end
    if  EstimOpt.Zinf == 1
        disp('Zero-inflated')
    end
    if isfield(EstimOpt,'Truncated') == 0 || EstimOpt.Truncated ==0
        disp('without truncation')
    elseif EstimOpt.Truncated == 1
        disp('with zero-truncation')    
    elseif EstimOpt.Truncated == 2
        disp('with zero-truncation and endogenous stratification (On-site sample)')
    end
    if isfield(EstimOpt,'Censored') == 0 || EstimOpt.Censored < 1
        disp('without censoring.')
    else
        disp(num2str(EstimOpt.Censored, 'with censoring at %2.0f.'))
    end
end

if isfield(EstimOpt,'Truncated') == 0 
    EstimOpt.Truncated = 0;
end

if isfield(EstimOpt,'Censored') == 0|| EstimOpt.Censored < 1
    EstimOpt.Censored = 0;
end

if EstimOpt.Zinf == 1
    if EstimOpt.NB == 0
        ModelTag = 'ZIP';
    elseif EstimOpt.NBtype == 1
        ModelTag = 'ZINB1';
    elseif EstimOpt.NBtype == 3
        ModelTag = 'ZINBP';
    else
        ModelTag = 'ZINB';
    end
elseif EstimOpt.NB == 0
    ModelTag = 'POISS';
elseif EstimOpt.NBtype == 1
    ModelTag = 'NB1';
elseif EstimOpt.NBtype == 3
    ModelTag = 'NBP';
else
    ModelTag = 'NB';
end

if ~ismember(EstimOpt.NB,[0 1])
    error('EstimOpt.NB must be 0 (Poisson) or 1 (Negative Binomial). For NB1 set EstimOpt.NB = 1 and EstimOpt.NBtype = 1')
end

if isfield(EstimOpt,'NBtype') == 0 || isempty(EstimOpt.NBtype)
    EstimOpt.NBtype = 2; % NB2 - variance lambda + alpha*lambda^2
elseif ~ismember(EstimOpt.NBtype,[1 2 3])
    error('EstimOpt.NBtype must be 1 (NB1, variance lambda*(1+alpha)), 2 (NB2, variance lambda + alpha*lambda^2) or 3 (NB-P, variance lambda + alpha*lambda^P with P estimated)')
end

% NB-P estimates the power P as a free parameter, placed right after the alpha block
EstimOpt.NVarNBP = double(EstimOpt.NB == 1 && EstimOpt.NBtype == 3);

if EstimOpt.Zinf == 1 && EstimOpt.Truncated > 0
    error('EstimOpt.Zinf cannot be combined with EstimOpt.Truncated - a zero-inflated model is not identified in a sample without zeros. Consider a hurdle specification instead.')
end

if isfield(EstimOpt,'HessEstFix') && EstimOpt.HessEstFix == 4 % analytical Hessian not derived for count data models yet
    if EstimOpt.Display == 1
        cprintf(rgb('DarkOrange'), 'WARNING: Analytical Hessian not available for count data models - using numerical Hessian instead (EstimOpt.HessEstFix = 3) \n')
    end
    EstimOpt.HessEstFix = 3;
end

% dependent variable
if any(~isfinite(INPUT.Y)) || any(INPUT.Y < 0) || any(INPUT.Y ~= floor(INPUT.Y))
    error('INPUT.Y must contain non-negative integers (counts)')
end

if EstimOpt.Censored > 0 && any(INPUT.Y > EstimOpt.Censored)
    if EstimOpt.Display == 1
        cprintf(rgb('DarkOrange'), ['WARNING: ',num2str(sum(INPUT.Y > EstimOpt.Censored)),' observation(s) with Y > EstimOpt.Censored - censoring them at ',num2str(EstimOpt.Censored),'. \n'])
    end
    INPUT.Y = min(INPUT.Y,EstimOpt.Censored);
end

if EstimOpt.Truncated > 0 && any(INPUT.Y == 0)
    error('Zero-truncated model requested (EstimOpt.Truncated > 0) but INPUT.Y contains %d zero(s)',sum(INPUT.Y == 0))
end

% explanatory variables
% Y was checked above; Xa was not checked at all before this, and a collinear
% design is the failure a travel cost package most needs to catch.  It does not
% stop the likelihood from converging - the fit is the same fit - but it splits
% the travel cost coefficient arbitrarily across the collinear columns, and
% consumer surplus is -1/beta read from that one slot.  See CDM_checkdesign.
Results.DesignCheck = CDM_checkdesign(INPUT,EstimOpt);

% weights - DataCleanCDM leaves them at the respondent level, the likelihood needs
% one per observation (cf. MNL.m). Never replace mis-sized weights silently.
NObs = size(INPUT.Y,1);
if isfield(INPUT,'W') == 0 || isempty(INPUT.W)
    INPUT.W = ones(NObs,1);
else
    INPUT.W = INPUT.W(:);
    if numel(INPUT.W) == NObs % already at the observation level
    elseif isfield(EstimOpt,'NCTMiss') && numel(INPUT.W) == numel(EstimOpt.NCTMiss)
        NCTMiss = EstimOpt.NCTMiss(:);
        Wtmp = zeros(NObs,1);
        Wtmp(1:NCTMiss(1)) = INPUT.W(1);
        for i = 2:numel(NCTMiss)
            Wtmp(sum(NCTMiss(1:i-1))+1:sum(NCTMiss(1:i))) = INPUT.W(i);
        end
        INPUT.W = Wtmp;
    else
        error('Incorrect size of the weights vector - expected %d (one per observation) or %d (one per respondent)',NObs,EstimOpt.NP)
    end
    INPUT.W = NObs*INPUT.W/sum(INPUT.W);
end
if isfield(EstimOpt,'NamesA') == 0 || isempty(EstimOpt.NamesA) || length(EstimOpt.NamesA) ~= EstimOpt.NVarA
    EstimOpt.NamesA = (1:EstimOpt.NVarA)';
    EstimOpt.NamesA = cellstr(num2str(EstimOpt.NamesA));
    if EstimOpt.Display == 1
        disp('Variable names (NamesA) unspecified, numbers used instead')
    end
elseif size(EstimOpt.NamesA,1) ~= EstimOpt.NVarA
    EstimOpt.NamesA = EstimOpt.NamesA';
end
% covariates of overdispersion - INPUT.Xnb gets a constant prepended, so the user
% supplies NVarNB-1 names and 'Cons' labels the intercept of the ln(alpha) equation
if EstimOpt.NVarNB > 1
    if isfield(EstimOpt,'NamesNB') == 0 || isempty(EstimOpt.NamesNB) || length(EstimOpt.NamesNB) ~= (EstimOpt.NVarNB-1)
        EstimOpt.NamesNB = (1:EstimOpt.NVarNB-1)';
        EstimOpt.NamesNB = cellstr(num2str(EstimOpt.NamesNB));
        if EstimOpt.Display == 1
            disp('Variable names (NamesNB) unspecified, numbers used instead')
        end
    elseif size(EstimOpt.NamesNB,1) ~= EstimOpt.NVarNB - 1
        EstimOpt.NamesNB = EstimOpt.NamesNB';
    end
    EstimOpt.NamesNB = [{'Cons'};EstimOpt.NamesNB];
else
    EstimOpt.NamesNB = {'ln alpha'};
end

if EstimOpt.NVarNBP == 1
    EstimOpt.NamesNB = [EstimOpt.NamesNB(:);{'P'}];
end

if EstimOpt.NVarZinf > 1 
    if isfield(EstimOpt,'NamesZinf') == 0 || isempty(EstimOpt.NamesZinf)|| length(EstimOpt.NamesZinf) ~= (EstimOpt.NVarZinf-1)
        EstimOpt.NamesZinf = (1:EstimOpt.NVarZinf-1)';
        EstimOpt.NamesZinf = cellstr(num2str(EstimOpt.NamesZinf));
    elseif size(EstimOpt.NamesZinf,1) ~= EstimOpt.NVarZinf - 1
        EstimOpt.NamesZinf = EstimOpt.NamesZinf';
    end
    EstimOpt.NamesZinf = [{'Cons'};EstimOpt.NamesZinf];
else
    EstimOpt.NamesZinf = {'Cons'};
end

%% starting values 

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + EstimOpt.NVarNB + EstimOpt.NVarNBP + EstimOpt.NVarZinf
    b0 = B_backup(:);
    if EstimOpt.Display == 1
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,ModelTag) && isfield(Results_old.(ModelTag),'b0') % starting values provided
    b0source = ModelTag;
    Results_old.(b0source).b0_old = Results_old.(b0source).b0;
    Results_old.(b0source) = rmfield(Results_old.(b0source),'b0');
    if length(Results_old.(b0source).b0_old) ~=  EstimOpt.NVarA + EstimOpt.NVarNB + EstimOpt.NVarNBP + EstimOpt.NVarZinf
       if EstimOpt.Display == 1
            disp('WARNING: Incorrect no. of starting values or model specification')
            disp(['Number of values in Results_old.',b0source,'.b0_old different from number of parameters for estimation'])
       end
       Results_old.(b0source) = rmfield(Results_old.(b0source),'b0_old');
    else
        b0 = Results_old.(b0source).b0_old(:);
    end
end


if  ~exist('b0','var')
    if EstimOpt.Display == 1
        disp('Using linear regression estimates as starting values')
    end
    b0 = INPUT.Xa\log(INPUT.Y+0.5); % OLS on log(Y+0.5), no Statistics Toolbox needed

    if EstimOpt.NVarNB > 0 % intercept plus any covariates of overdispersion
        b0 = [b0; zeros(EstimOpt.NVarNB,1)];
    end

    if EstimOpt.NVarNBP == 1
        b0 = [b0; 2]; % start NB-P at P = 2, i.e. at NB2
    end

    if EstimOpt.Zinf == 1
        b0 = [b0; zeros(EstimOpt.NVarZinf,1)];
    end
end

%% Optimization Options
if isfield(EstimOpt,'BActive')
    EstimOpt.BActive = EstimOpt.BActive(:)';
end

if EstimOpt.ConstVarActive == 1
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        error ('Are there any constraints on model parameters (EstimOpt.ConstVarActive)? Constraints not provided (EstimOpt.BActive).')
    elseif length(b0) ~= length(EstimOpt.BActive)
        error('Check no. of constraints, length(b0) ~= length(EstimOpt.BActive)')
    end
    if EstimOpt.Display == 1
        disp(['Initial values: ' mat2str(b0',2)])
        disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
    end
else
    if ~isfield(EstimOpt,'BActive') || isempty(EstimOpt.BActive) || sum(EstimOpt.BActive == 0) == 0
        EstimOpt.BActive = ones(1,length(b0));
        if EstimOpt.Display == 1
            disp(['Initial values: ' mat2str(b0',2)])
        end
    else
        if length(b0) ~= length(EstimOpt.BActive)
            error('Check no. of constraints, length(b0) ~= length(EstimOpt.BActive)')
        else
            if EstimOpt.Display == 1
                disp(['Initial values: ' mat2str(b0',2)])
                disp(['Parameters with zeros are constrained to their initial values: ' mat2str(EstimOpt.BActive')])
            end
        end
    end
end

if ((isfield(EstimOpt, 'ConstVarActive') == 1 && EstimOpt.ConstVarActive == 1) || sum(EstimOpt.BActive == 0) > 0) && ~isequal(OptimOpt.GradObj,'on')
    if EstimOpt.Display == 1
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied gradient on - otherwise parameters'' constraints will be ignored - switch to constrained optimization instead (EstimOpt.ConstVarActive = 1) \n')
    end
    OptimOpt.GradObj = 'on';
end

if (isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0) && isequal(OptimOpt.Algorithm,'quasi-newton') && isequal(OptimOpt.Hessian,'user-supplied')
    if EstimOpt.Display == 1
        cprintf(rgb('DarkOrange'), 'WARNING: Setting user-supplied Hessian off - quasi-newton algorithm does not use it anyway \n')
    end
    OptimOpt.Hessian = 'off';
end

if EstimOpt.Display == 1
    fprintf('\n')
    cprintf('Opmization algorithm: '); cprintf('*Black',[OptimOpt.Algorithm '\n'])
    if strcmp(OptimOpt.GradObj,'on')
        if EstimOpt.NumGrad == 0
            cprintf('Gradient: '); cprintf('*Black','user-supplied, analytical \n')
        else
            cprintf('Gradient: '); cprintf('*Black',['user-supplied, numerical, ' OptimOpt.FinDiffType '\n'])
        end
    else
        cprintf('Gradient: '); cprintf('*Black',['built-in, ' OptimOpt.FinDiffType '\n'])
    end
%from Wiktors code, possibility to switch how is hessian included
    if isequal(OptimOpt.Algorithm,'quasi-newton') 
        cprintf('Hessian: '); cprintf('*Black','off, ')
        switch EstimOpt.HessEstFix
            case 0
                cprintf('*Black','retained from optimization \n')
            case 1
                cprintf('*Black','ex-post calculated using BHHH \n')
            case 2
                cprintf('*Black','ex-post calculated using high-precision BHHH \n')
            case 3
                cprintf('*Black','ex-post calculated numerically \n')
            case 4
                cprintf('*Black','ex-post calculated analytically \n')
        end
    else
        if strcmp(OptimOpt.Hessian,'user-supplied')
            cprintf('Hessian: '); cprintf('*Black','user-supplied, BHHH, ')
        else
            cprintf('Hessian: '); cprintf('*Black',['built-in, ' OptimOpt.HessUpdate ', '])
        end
        switch EstimOpt.HessEstFix
            case 0
                cprintf('*Black','retained from optimization \n')
            case 1
                cprintf('*Black','ex-post calculated using BHHH \n')
            case 2
                cprintf('*Black','ex-post calculated using high-precision BHHH \n')
            case 3
                cprintf('*Black','ex-post calculated numerically \n')
            case 4
                cprintf('*Black','ex-post calculated analytically \n')
        end
    end
    fprintf('\n')
end

%% estimation
%LLfun = @(B) LL_CDM_MATlike(INPUT.Y, INPUT.Xa,INPUT.Xnb,INPUT.W, EstimOpt,OptimOpt,B);
LLfun = @(B) LL_CDM_MATlike(INPUT.Y, INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,INPUT.W, EstimOpt,OptimOpt,B);


if EstimOpt.ConstVarActive == 0 
    
    if EstimOpt.HessEstFix == 0 %optimization hessian used
        [Results.bhat, LL, Results.exitf, Results.output, Results.g, Results.hess] = fminunc(LLfun, b0, OptimOpt);
    else
        [Results.bhat, LL, Results.exitf, Results.output, Results.g] = fminunc(LLfun, b0, OptimOpt);
    end  
    
elseif EstimOpt.ConstVarActive == 1 % equality constraints
        
    EstimOpt.CONS1 = diag(1 - EstimOpt.BActive);
    EstimOpt.CONS1(sum(EstimOpt.CONS1,1)==0,:)=[];
    EstimOpt.CONS2 = zeros(size(EstimOpt.CONS1,1),1);

    if EstimOpt.HessEstFix == 0
        [Results.bhat, LL, Results.exitf, Results.output, Results.lambda, Results.g, Results.hess] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    else
        [Results.bhat, LL, Results.exitf, Results.output, Results.lambda, Results.g] = fmincon(LLfun,b0,[],[],EstimOpt.CONS1,EstimOpt.CONS2,[],[],[],OptimOpt);
    end

end

%% Saving results

Results.LL = -LL;
Results.b0_old = b0;

% LL_CDM returns the negative log-likelihood per observation; store the weighted
% log-likelihood, which is what the LL reported above is the sum of.
Results.LLdetailed = -INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);

if EstimOpt.HessEstFix == 1
    Results.jacobian = numdiff(@(B) INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),-Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),1), Results.bhat);
end


if sum(EstimOpt.BActive == 0) > 0
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.jacobian = Results.jacobian(:, EstimOpt.BActive == 1);
        Results.hess = Results.jacobian'*Results.jacobian;
    elseif EstimOpt.HessEstFix == 0 || EstimOpt.HessEstFix == 3
        Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
    end
    [Results.ihess, Results.HessDiagnostics] = hessianInverse(Results.hess,ModelTag);
    Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
    Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
    
else
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.hess = Results.jacobian'*Results.jacobian;
    end
    [Results.ihess, Results.HessDiagnostics] = hessianInverse(Results.hess,ModelTag);
end


if EstimOpt.RobustStd == 1
    if EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),-Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
    end
    RobJacob = zeros(EstimOpt.NP, size(Results.jacobian,2));
    RobJacob(1,:) = sum(Results.jacobian(1:EstimOpt.NCTMiss(1),:),1);
    for i = 2:EstimOpt.NP
        RobJacob(i,:) = sum(Results.jacobian(sum(EstimOpt.NCTMiss(1:(i-1)))+1:sum(EstimOpt.NCTMiss(1:i)),:),1);
    end
    RobustHess = RobJacob'*RobJacob;
    Results.ihess = Results.ihess*RobustHess*Results.ihess;
end
Results.std = sqrt(diag(Results.ihess));


if sum(EstimOpt.BActive == 0) > 0
    Results.std(EstimOpt.BActive == 0) = NaN;
end
Results.std(imag(Results.std) ~= 0) = NaN;
Results.R = [Results.bhat , Results.std , pv(Results.bhat , Results.std)];

EstimOpt.Params = length(b0);
if isfield(EstimOpt,'BActive')
    EstimOpt.Params = EstimOpt.Params - sum(EstimOpt.BActive == 0);
end

if isfield(Results_old,'CDM0') && isfield(Results_old.CDM0,'LL')
    LL0 = Results_old.CDM0.LL;
elseif isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL') % pre-2026 field name
    LL0 = Results_old.POISS0.LL;
else
    LL0 = [];
end

if ~isempty(LL0)
    % slot 4 is Ben-Akiva-Lerman's pseudo-R2, which is not defined for count models
    Results.stats = [Results.LL; LL0; 1-Results.LL/LL0; NaN; ((2*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs; ((log(EstimOpt.NObs)*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs; EstimOpt.NObs; EstimOpt.NP; EstimOpt.Params];
else
    Results.stats = [];
end

Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;

Results.DetailsA = zeros(EstimOpt.NVarA,4);
Results.DetailsA(:,[1 3 4]) = [Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))] ;

if EstimOpt.NB == 1
    Results.DetailsNB = zeros(EstimOpt.NVarNB+EstimOpt.NVarNBP,4);
    Results.DetailsNB(:,[1 3 4]) = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP))] ;
end

if EstimOpt.Zinf == 1
    Results.DetailsZinf = zeros(EstimOpt.NVarZinf,4);
    Results.DetailsZinf(:,[1 3 4]) = [Results.bhat(EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+EstimOpt.NVarZinf),Results.std(EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+EstimOpt.NVarZinf),pv(Results.bhat(EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+EstimOpt.NVarZinf),Results.std(EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarNBP+EstimOpt.NVarZinf))] ;
end

if TC_par_index > 0 && TC_par_index <= EstimOpt.NVarA
    Results.DetailsCS = zeros(1,4);
    Results.DetailsCS(:,[1 3 4]) = [-1/Results.DetailsA(TC_par_index,1), sqrt(Results.ihess(TC_par_index,TC_par_index)*((-1./Results.DetailsA(TC_par_index,1))^4)), pv(-1/Results.DetailsA(TC_par_index,1), sqrt(Results.ihess(TC_par_index,TC_par_index)*((-1./Results.DetailsA(TC_par_index,1))^4)))];
end


%% Template filling


Template1 = {'DetailsA'};
Template2 = {'DetailsA'};
Names.DetailsA = EstimOpt.NamesA;
Heads.DetailsA = {'';'tb'};
ST = {'DetailsA'};

if EstimOpt.NVarNB > 0
    Template1 = [Template1;{'DetailsNB'}];
    Template2 = [Template2;{'DetailsNB'}];
    Names.DetailsNB = EstimOpt.NamesNB;
    Heads.DetailsNB(1,1:2) = {'','Overdispersion parameters'};
    Heads.DetailsNB(end+1,1:2) = {'tb','tb'};
    ST = [ST,{'DetailsNB'}];
end

if EstimOpt.NVarZinf > 0
    Template1 = [Template1;{'DetailsZinf'}];
    Template2 = [Template2;{'DetailsZinf'}];
    Names.DetailsZinf = EstimOpt.NamesZinf;
    Heads.DetailsZinf(1,1:2) = {'','Zero inflation'};
    Heads.DetailsZinf(end+1,1:2) = {'tb','tb'};
    ST = [ST,{'DetailsZinf'}];
end

if isfield(Results,'DetailsCS')
    Template1 = [Template1;{'DetailsCS'}];
    Template2 = [Template2;{'DetailsCS'}];
    Names.DetailsCS = cellstr(strcat('-1/(', EstimOpt.NamesA(TC_par_index),')'));
    Heads.DetailsCS(1,1:2) = {'','Consumer Surplus'};
    Heads.DetailsCS(end+1,1:2) = {'tb','tb'};
    ST = [ST,{'DetailsCS'}];
end

%% Header

Head = cell(1,2);

if EstimOpt.NB == 0
    if EstimOpt.Zinf == 1
        Head(1,1) = {'Zero-inflated Poisson'};
    else
        Head(1,1) = {'Poisson'};
    end
else
    if EstimOpt.NBtype == 1
        nbName = 'Negative Binomial (NB1)';
    elseif EstimOpt.NBtype == 3
        nbName = 'Negative Binomial (NB-P)';
    else
        nbName = 'Negative Binomial';
    end
    if EstimOpt.Zinf == 1
        Head(1,1) = {['Zero-inflated ',nbName]};
    else
        Head(1,1) = {nbName};
    end
end
if isfield(EstimOpt,'Censored') == 0 || EstimOpt.Censored < 1
    if EstimOpt.Truncated == 1
        Head(1,2) = {'with zero-truncation'};
    elseif EstimOpt.Truncated == 2
        Head(1,2) = {'with zero-truncation and endogenous stratification'};
    else
        Head(1,2) = {''};
    end
else 
    if EstimOpt.Truncated == 1
        Head(1,2) = {num2str(EstimOpt.Censored, 'with zero-truncation, censored at %2.0f')};
    elseif EstimOpt.Truncated == 2
        Head(1,2) = {num2str(EstimOpt.Censored,'with zero-truncation and endogenous stratification, censored at %2.0f')};
    else
        Head(1,2) = {num2str(EstimOpt.Censored,'censored at %2.0f')};
    end
end


%% Footer

% Rows 1-16 follow the layout genOutput expects (as in MNL.m); the descriptive
% statistics table is appended below it, in the extra columns.
Tail = cell(16,2);
Tail(2,1) = {'Model diagnostics'};
Tail(3:16,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));strcat('Ben-Akiva-Lerman''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';' ';'Estimation method';'Optimization method';'Gradient';'Hessian'};

if ~isempty(Results.stats)
    Tail(3:11,2) = num2cell(Results.stats);
    Tail(6,2) = {[]}; % Ben-Akiva-Lerman's pseudo-R2 is not defined for count models
end


if any(INPUT.W ~= 1)
    Tail(13,2) = {'weighted maximum likelihood'};
else
    Tail(13,2) = {'maximum likelihood'};
end

Tail(14,2) = {OptimOpt.Algorithm;};

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        Tail(15,2) = {'user-supplied, analytical'};
    else
        Tail(15,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
    end
else
    Tail(15,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
end

if isequal(OptimOpt.Algorithm,'quasi-newton')
    outHessian='off, ';
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian,'retained from optimization'];
        case 1
            outHessian = [outHessian,'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian,'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian,'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian,'ex-post calculated analytically'];
    end
else
    if strcmp(OptimOpt.Hessian,'user-supplied')
        if EstimOpt.ApproxHess == 1
            outHessian = 'user-supplied, BHHH, ';
        else
            outHessian = 'user-supplied, analytical, ';
        end
    else
        outHessian = ['built-in, ',num2str(OptimOpt.HessUpdate),', '];
    end
    switch EstimOpt.HessEstFix
        case 0
            outHessian = [outHessian,'retained from optimization'];
        case 1
            outHessian = [outHessian,'ex-post calculated using BHHH'];
        case 2
            outHessian = [outHessian,'ex-post calculated using high-precision BHHH'];
        case 3
            outHessian = [outHessian,'ex-post calculated numerically'];
        case 4
            outHessian = [outHessian,'ex-post calculated analytically'];
    end
end

Tail(16,2) = {outHessian};

if TC_par_index > 0 && TC_par_index <= EstimOpt.NVarA % widens Tail to 3 columns
    Tail(18,1) = {'Descriptive statistics'};
    Tail(19,1:3) = {'var','Y',char(EstimOpt.NamesA(TC_par_index))};
    desc_y = dstats(INPUT.Y);
    desc_tc = dstats(INPUT.Xa(:,TC_par_index));
    Tail(20:26,1:3) = [desc_y(1:7,:),desc_tc(1:7,2)];
end

%%  Print to screen and .xls

if EstimOpt.Display ~= 0
    Results.Dist = -ones(EstimOpt.NVarA,1);
    if ~isfield(EstimOpt,'OutputFilePrefix') || isempty(EstimOpt.OutputFilePrefix)
        % Head{1,1} is a sentence ('Zero-inflated Negative Binomial'); the saved file
        % gets a short tag instead, in line with MNL_results_*.xls in DCE.
        if isfield(EstimOpt,'ProjectName') && ~isempty(EstimOpt.ProjectName)
            EstimOpt.OutputFilePrefix = [ModelTag,'_results_',char(string(EstimOpt.ProjectName))];
        else
            EstimOpt.OutputFilePrefix = [ModelTag,'_results'];
        end
    end
    if EstimOpt.NB == 1
        % H0: ln(alpha) = 0 is not the Poisson boundary, so the intercept of the
        % overdispersion equation carries no stars. Its covariates are ordinary tests.
        if EstimOpt.NVarNB == 1
            EstimOpt.NoSigRows.DetailsNB = Inf; % the block is the intercept alone
        else
            EstimOpt.NoSigRows.DetailsNB = 1;
        end
    end
    Results.R_out = genOutput(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);
end
end


function X = addConstant(X,name,Display)
% addConstant prepends a column of ones, dropping one the user already supplied.
% CDM owns the constant for Xnb and Xzinf, so its label ('Cons') always lines up
% with the first column and NamesNB/NamesZinf list the covariates only.
isConst = all(X == X(1,:),1);
if any(isConst)
    if Display == 1
        cprintf(rgb('DarkOrange'), ['WARNING: INPUT.',name,' contains a constant column - CDM adds its own, so the supplied one is dropped \n'])
    end
    X(:,isConst) = [];
end
X = [ones(size(X,1),1), X];
end


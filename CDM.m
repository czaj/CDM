function Results = CDM(INPUT,Results_old,EstimOpt,OptimOpt)
% POISS creates a simple count model
%   That is either Poisson, Negative binomial 
%
% Syntax:   POISS(INPUT,EstimOpt,OptimOpt)
%           POISS(INPUT,Results_old,EstimOpt,OptimOpt)
%
% Inputs:
%    INPUT - clean, updated INPUT data from DataCleanTCM
%    EstimOpt - Estimation Options (check below)
%    OptimOpt - Optimizer Options define how algorithms converges to the final result. 
%    They are set by default based on provided EstimOpt in DataCleanTCM, however, they are subject to change.
%    Results_old - here one can provide old results to use as starting
%    values
%
% EstimOpt Options:
% Set them by e.g. Estimopt.DataFile = 'Project'
%
% General basics:
% 	DataFile  path/name of the .mat data file
% 	Display  1; shows output, set to 0 to hide it 
% 	ProjectName  Name of the project/model
% 	WTP_space  set to 1 for estimation in WTP space. If missing or set to 0, MNL uses Preference Space
% 	NCT - Number of choice tasks per person 
% 	NAlt - Number of alternatives
% 	NP  Number of respondents
% 
% 
% Variables options:
% 	NamesA  Names of variables in list e.g. {'-Opt out';-Cost (EUR)'}
% 	NamesM  Names of variables of means of random parameters
% 	NamesS  Names of variables of Scale
% 	NLTVariables  vector specifying which attributes are to subject to non-linear transformations
% 	NLTType  Transformation for non-linear variables. By default it is set to Box-Cox transformation (1), set to 2 in order to use Yeo-Johnson transformation
% 
% 
% Parameters options:
% 	ExpB = vector of 0; for each parameter set it to 1 to use ExpB, otherwise 0
% 	BActive = vector of 0; for each parameter set it to 1 to constrain model parameters to their initial values
% 	ConstVarActive = 0; set to 1 to constrain model parameters to its initial values 
% 
% 
% Modelling options from DataCleanDCE:
% 	ApproxHess = 1; for user supplied hessians, 1 for BHHH, 0 for analytical
% 	RobustStd = 0; by default not using robust standard errors, set to 1 to use them
% 	NumGrad = 0; uses analytical gradient in calculations, set to 1 for numerical gradient
% 	HessEstFix = 0; Options: 
% o	0 - use optimization Hessian, 
% o	1 - use jacobian-based (BHHH) Hessian, 
% o	2 - use high-precision jacobian-based (BHHH) Hessian,
% o	3 - use numerical Hessian, 
% o	4 - use analytical Hessian
%
% Example: 
%    Results.POISS = POISS(INPUT,Results,EstimOpt,OptimOpt);


%TC_par_index - which parameter to use for CS calculation, default = 2
%clear
%clc
global B_backup;

%% 
tic 

Results.bhat = [];
Results.R = [];
%Results.R_out = {};
Results.stats = [];

%% check input:

if nargin<3
    error('Too few input arguments for POISS(INPUT,EstimOpt,OptimOpt)')
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

if isfield(EstimOpt, 'TC_par_index') == 0 || EstimOpt.TC_par_index ==0
    EstimOpt.TC_par_index = 2;
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
    INPUT.Xnb = [ones(size(INPUT.Y,1),1), INPUT.Xnb];
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
    INPUT.Xzinf = INPUT.Xzinf;
    disp('test done');
end

EstimOpt.NVarZinf = size(INPUT.Xzinf,2);

if EstimOpt.Display == 1
    disp(' ');
    disp('__________________________________________________________________________________________________________________');
    disp(' ');
    if EstimOpt.NB == 0
        disp('Estimating Poisson regression model ...')
    elseif EstimOpt.NB == 1
        disp('Estimating Negative Binomial regression model ...')
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

% weights
if isfield(INPUT,'W') == 0 || length(INPUT.W) ~= length(INPUT.Y)
    INPUT.W = ones(length(INPUT.Y),1);
else
    INPUT.W = INPUT.W(:);
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
%if add NB>1
if EstimOpt.NVarNB > 1 
    % if isfield(EstimOpt,'NamesNB') == 0 || isempty(EstimOpt.NamesNB)|| length(EstimOpt.NamesNB) ~= (EstimOpt.NVarNB-1)
    %     EstimOpt.NamesNB = (1:EstimOpt.NVarNB-1)';
    %     EstimOpt.NamesNB = cellstr(num2str(EstimOpt.NamesNB));
    % elseif size(EstimOpt.NamesNB,1) ~= EstimOpt.NVarNB - 1
    %     EstimOpt.NamesNB = EstimOpt.NamesNB';
    % end
    % EstimOpt.NamesNB = [{'Cons'};EstimOpt.NamesNB];
else
    EstimOpt.NamesNB = {'ln alpha'};
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

if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + EstimOpt.NVarNB + EstimOpt.NVarZinf
    b0 = B_backup(:);
    if EstimOpt.Display == 1
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,'POISS') && isfield(Results_old.POISS,'b0') % starting values provided

    Results_old.POISS.b0_old = Results_old.POISS.b0;
    Results_old.POISS = rmfield(Results_old.POISS,'b0');
    if length(Results_old.POISS.b0_old) ~=  EstimOpt.NVarA + EstimOpt.NVarNB + EstimOpt.NVarZinf
       if EstimOpt.Display == 1
            disp('WARNING: Incorrect no. of starting values or model specification')
            disp('Number of values in Results_old.POISS.b0_old different from number of parameters for estimation')
       end
       Results_old.POISS = rmfield(Results_old.POISS,'b0_old');
    else
        b0 = Results_old.POISS.b0_old(:);
    end
end


if  ~exist('b0','var')
    if EstimOpt.Display == 1
        disp('Using linear regression estimates as starting values')
    end
    b0 = regress(log(INPUT.Y+0.1), INPUT.Xa);

    if EstimOpt.NVarNB == 1
        b0 = [b0; zeros(EstimOpt.NVarNB,1)];
    end

    %NOWY KOD
    %jakie miejsca startowe?

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
        [Results.bhat, LL, Results.exitf, Results.output, Results.g, Results.hess] = fminunc(LLfun, b0, OptimOpt); %tu się dzieje coś dziwnego z b0, dlaczego w kolejnych iteracjach std. err maleje?
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

if EstimOpt.HessEstFix == 1
	f = LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf, EstimOpt,B),f,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),1), Results.bhat);
end


if sum(EstimOpt.BActive == 0) > 0
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.jacobian = Results.jacobian(:, EstimOpt.BActive == 1);
        Results.hess = Results.jacobian'*Results.jacobian;
    elseif EstimOpt.HessEstFix == 0 || EstimOpt.HessEstFix == 3
        Results.hess = Results.hess(EstimOpt.BActive == 1,EstimOpt.BActive == 1);
    end
    Results.ihess = inv(Results.hess);
    Results.ihess = direcXpnd(Results.ihess,EstimOpt.BActive);
    Results.ihess = direcXpnd(Results.ihess',EstimOpt.BActive);
    
else
    if EstimOpt.HessEstFix == 1 || EstimOpt.HessEstFix == 2
        Results.hess = Results.jacobian'*Results.jacobian;
    end
    Results.ihess = inv(Results.hess);
end


if EstimOpt.RobustStd == 1
    if EstimOpt.NumGrad == 0
        [~, Results.jacobian] = LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_CDM(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
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

if isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL')
    Results.stats = [Results.LL; Results_old.POISS0.LL;  1-Results.LL/Results_old.POISS0.LL; NaN; ((2*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs; ((log(EstimOpt.NObs)*EstimOpt.Params-2*Results.LL))/EstimOpt.NObs ;EstimOpt.NObs; EstimOpt.NP; EstimOpt.Params];

else
    Results.stats = NaN(9,1);
end

Results.INPUT = INPUT;
Results.EstimOpt = EstimOpt;
Results.OptimOpt = OptimOpt;

Results.DetailsA = zeros(EstimOpt.NVarA,4);
Results.DetailsA(:,[1 3 4]) = [Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA),pv(Results.bhat(1:EstimOpt.NVarA),Results.std(1:EstimOpt.NVarA))] ;

if EstimOpt.NB == 1
    Results.DetailsNB = zeros(EstimOpt.NVarNB,4);
    Results.DetailsNB(:,[1 3 4]) = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB))] ;
elseif EstimOpt.NB == 2
    Results.DetailsNB = zeros(EstimOpt.NVarNB+1,4);
    Results.DetailsNB(:,[1 3 4]) = [Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+1),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+1),pv(Results.bhat(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+1),Results.std(EstimOpt.NVarA+1:EstimOpt.NVarA+EstimOpt.NVarNB+1))] ;   
end

if EstimOpt.Zinf == 1
    Results.DetailsZinf = zeros(EstimOpt.NVarZinf,4);
    Results.DetailsZinf(:,[1 3 4]) = [Results.bhat(EstimOpt.NVarA+EstimOpt.NVarNB+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarZinf),Results.std(EstimOpt.NVarA+EstimOpt.NVarNB+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarZinf),pv(Results.bhat(EstimOpt.NVarA+EstimOpt.NVarNB+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarZinf),Results.std(EstimOpt.NVarA+EstimOpt.NVarNB+1:EstimOpt.NVarA+EstimOpt.NVarNB+EstimOpt.NVarZinf))] ;
end

if isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL')
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

if isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL')
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
    if EstimOpt.Zinf == 1
        Head(1,1) = {'Zero-inflated Negative Binomial'};
    else
        Head(1,1) = {'Negative Binomial'};
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

Tail = cell(26,4);
Tail(2,1) = {'Model diagnostics'};
Tail(3:15,1) = {'LL at convergence';'LL at constant(s) only';strcat('McFadden''s pseudo-R',char(178));'AIC/n';'BIC/n';'n (observations)';'r (respondents)';'k (parameters)';' ';'Estimation method';'Optimization method';'Gradient';'Hessian'};

if isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL')
    Tail(3:10,2) = num2cell(Results.stats([1:3,5:9]));
end


if any(INPUT.W ~= 1)
    Tail(12,2) = {'weighted maximum likelihood'};
else
    Tail(12,2) = {'maximum likelihood'};
end

Tail(13,2) = {OptimOpt.Algorithm;};

if strcmp(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        Tail(14,2) = {'user-supplied, analytical'};
    else
        Tail(14,2) = {['user-supplied, numerical ',num2str(OptimOpt.FinDiffType)]};
    end
else
    Tail(14,2) = {['built-in, ',num2str(OptimOpt.FinDiffType)]};
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

Tail(15,2) = {outHessian};

Tail(17,1) = {'Descriptive Statistics'};

if isfield(Results_old,'POISS0') && isfield(Results_old.POISS0,'LL')
    Tail(18,1:3) = {'var','Y',char(EstimOpt.NamesA(TC_par_index))};
    desc_y = dstats(INPUT.Y);
    desc_tc = dstats(INPUT.Xa(:,TC_par_index));
    desc = [desc_y(1:7,:),desc_tc(1:7,2)];
    Tail(20:26,1:3) = desc;
end

%%  Print to screen and .xls

if EstimOpt.Display ~= 0
    Results.Dist = -ones(EstimOpt.NVarA,1);
    Results.R_out = genOutput_CDM(EstimOpt,Results,Head,Tail,Names,Template1,Template2,Heads,ST);

end
end
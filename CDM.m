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

if isfield(EstimOpt, 'NB') == 0 
    EstimOpt.NB = 0;
end

if isfield(EstimOpt, 'Zinf') == 0 
    EstimOpt.Zinf = 0;
end

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
        disp('Zero inflated')
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
%!!!!!!!!!
if EstimOpt.NVarNB > 1 
    % if isfield(EstimOpt,'NamesNB') == 0 || isempty(EstimOpt.NamesNB)|| length(EstimOpt.NamesNB) ~= (EstimOpt.NVarNB-1)
    %     EstimOpt.NamesNB = (1:EstimOpt.NVarNB-1)';
    %     EstimOpt.NamesNB = cellstr(num2str(EstimOpt.NamesNB));
    % elseif size(EstimOpt.NamesNB,1) ~= EstimOpt.NVarNB - 1
    %     EstimOpt.NamesNB = EstimOpt.NamesNB';
    % end
    % EstimOpt.NamesNB = [{'Cons'};EstimOpt.NamesNB];
else
    EstimOpt.NamesNB = {'Cons'};
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

%if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + EstimOpt.NVarNB
if exist('B_backup','var') && ~isempty(B_backup) && size(B_backup,1) == EstimOpt.NVarA + EstimOpt.NVarNB + EstimOpt.NVarZinf
    b0 = B_backup(:);
    if EstimOpt.Display == 1
        disp('Using the starting values from Backup')
    end
elseif isfield(Results_old,'POISS') && isfield(Results_old.POISS,'b0') % starting values provided

    Results_old.POISS.b0_old = Results_old.POISS.b0;
    Results_old.POISS = rmfield(Results_old.POISS,'b0');
    %if length(Results_old.POISS.b0_old) ~=  EstimOpt.NVarA + EstimOpt.NVarNB
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

if EstimOpt.HessEstFix == 1
	f = LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
    Results.jacobian = numdiff(@(B) LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf, EstimOpt,B),f,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
elseif EstimOpt.HessEstFix == 2
    Results.jacobian = jacobianest(@(B) LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),Results.bhat);
elseif EstimOpt.HessEstFix == 3
    Results.hess = hessian(@(B) sum(LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),1), Results.bhat);
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
        [~, Results.jacobian] = LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
        Results.jacobian = Results.jacobian.*INPUT.W(:, ones(1,size(Results.jacobian,2)));
    else
        Results.LLdetailed = LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,Results.bhat);
        Results.jacobian = numdiff(@(B) INPUT.W.*LL_DCE(INPUT.Y,INPUT.Xa,INPUT.Xnb,INPUT.Xzinf,EstimOpt,B),Results.LLdetailed,Results.bhat,isequal(OptimOpt.FinDiffType, 'central'),EstimOpt.BActive);
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

if EstimOpt.Display == 1
    disp(' ');
    disp('Lambda parameters');
    disp(['var.', blanks(size(char(EstimOpt.NamesA),2)-2) ,'coef.      st.err.  p-value'])
    disp([char(EstimOpt.NamesA) ,blanks(EstimOpt.NVarA)', num2str(Results.DetailsA(:,1),'%8.4f') star_sig(Results.DetailsA(:,4)) num2str(Results.DetailsA(:,3:4),'%8.4f %8.4f')])

    if EstimOpt.NB  > 0 
        disp(' ');
        disp('Overdispersion parameters');
        disp(['var.', blanks(size(char(EstimOpt.NamesNB),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesNB) ,blanks(EstimOpt.NVarNB + (EstimOpt.NB == 2))', num2str(Results.DetailsNB(:,1),'%8.4f') star_sig(Results.DetailsNB(:,4)) num2str(Results.DetailsNB(:,3:4),'%8.4f %8.4f')])
    end

    if EstimOpt.Zinf  == 1 
        disp(' ');
        disp('Zero inflation ');
        disp(['var.', blanks(size(char(EstimOpt.NamesZinf),2)-2) ,'coef.      st.err.  p-value'])
        disp([char(EstimOpt.NamesZinf) ,blanks(EstimOpt.NVarZinf)', num2str(Results.DetailsZinf(:,1),'%8.4f') star_sig(Results.DetailsZinf(:,4)) num2str(Results.DetailsZinf(:,3:4),'%8.4f %8.4f')])
    end
end
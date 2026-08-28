function [INPUT, Results, EstimOpt, OptimOpt] = DataCleanCDM(INPUT,EstimOpt)
% DataCleanCDM prepares the data and options for CDM.
%
% It converts INPUT to double, removes observations flagged in INPUT.MissingInd,
% puts the weights on a per-respondent basis and scales them to unit mean, fills
% in the estimation and optimizer defaults, and estimates the constants-only
% model whose log-likelihood is used for McFadden's pseudo-R2.
%
% Syntax:   [INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt)
%
% Required EstimOpt fields: NP (respondents) and NCT (observations per
% respondent). See the help of CDM for everything else.
%
% See also CDM, LL_CDM, CDM_demo
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#


% global TolB
% save tmp1
% return

inputnames = fieldnames(INPUT);
for i=1:length(inputnames)
    INPUT.(inputnames{(i)}) = double(INPUT.(inputnames{(i)}));
end

if isfield(INPUT,'MissingInd') == 0 || isempty(INPUT.MissingInd)
    INPUT.MissingInd = zeros(size(INPUT.Y,1),1);
end
INPUT.MissingInd = double(INPUT.MissingInd(:) ~= 0);

EstimOpt.MissingCT = [];
NP_in = EstimOpt.NP; % no. of respondents before dropping the incomplete ones

% Weights are kept at the respondent level throughout DataCleanCDM (as in DataCleanDCE)
% and expanded to the observation level by the model function.
if isfield(INPUT,'W') && ~isempty(INPUT.W)
    INPUT.W = INPUT.W(:);
    if numel(INPUT.W) == EstimOpt.NCT*NP_in
        INPUT.W = INPUT.W(1:EstimOpt.NCT:end); % one weight per respondent
    elseif numel(INPUT.W) ~= NP_in
        error('Incorrect size of the weights vector - provide one weight per respondent (EstimOpt.NP) or per observation (EstimOpt.NCT*EstimOpt.NP)')
    end
end

if any(INPUT.MissingInd)
    MissingCT = reshape(INPUT.MissingInd,EstimOpt.NCT,NP_in) == 1; % missing observations
    MissingP = sum(MissingCT,1) == EstimOpt.NCT; % respondents with all NCT missing

    cprintf('text',['The dataset contains ',num2str(sum(INPUT.MissingInd)),' observation(s) with missing data (out of the total of ',num2str(numel(INPUT.MissingInd)),') - removing them.\n'])

    if any(MissingP) % respondents with 0 observations left
        cprintf(rgb('DarkOrange'),['WARNING: Dataset includes ',num2str(sum(MissingP)),' respondents with 0 completed observations. Adjusting NP from ',num2str(NP_in),' to ',num2str(NP_in-sum(MissingP)),'.\n'])
        EstimOpt.NP = NP_in - sum(MissingP);
        if isfield(INPUT,'W') && ~isempty(INPUT.W)
            INPUT.W = INPUT.W(~MissingP);
        end
    end

    EstimOpt.NCTMiss = EstimOpt.NCT - sum(MissingCT(:,~MissingP),1)';

    % Drop the incomplete observations from every observation-level field. Leaving
    % them in as NaN would make the log-likelihood NaN for the whole sample.
    keepRows = INPUT.MissingInd == 0;
    INPUT_fields = fieldnames(INPUT);
    for i = 1:numel(INPUT_fields)
        if strcmp(INPUT_fields{i},'W') % weights stay at the respondent level
            continue
        end
        tmp = INPUT.(INPUT_fields{i});
        if size(tmp,1) == EstimOpt.NCT*NP_in
            INPUT.(INPUT_fields{i}) = tmp(keepRows,:);
        end
    end
else
    EstimOpt.NCTMiss = EstimOpt.NCT*ones(EstimOpt.NP,1);
end

EstimOpt.NObs = sum(EstimOpt.NCTMiss);

if isfield(INPUT,'W') && ~isempty(INPUT.W)
    if abs(sum(INPUT.W) - EstimOpt.NP) > 1e-10*EstimOpt.NP
        cprintf(rgb('DarkOrange'), ['WARNING: Scaling weights for unit mean. \n'])
        INPUT.W = EstimOpt.NP*INPUT.W/sum(INPUT.W);
    end
else
    INPUT.W = ones(EstimOpt.NP,1);
end

if isfield(EstimOpt,'RobustStd') == 0
    EstimOpt.RobustStd = 0; % do not use robust standard errors
end

EstimOpt.NVarA = size(INPUT.Xa,2); % Number of attributes

%EstimOpt.NVarZinf = size(INPUT.Xzinf,2);

if isfield(EstimOpt,'ConstVarActive') == 0
    EstimOpt.ConstVarActive = 0;
end
if isfield(EstimOpt,'Display') == 0
    EstimOpt.Display = 1;
end

if isfield(EstimOpt,'NumGrad') == 0 || (EstimOpt.NumGrad ~= 0 && EstimOpt.NumGrad ~= 1)
    EstimOpt.NumGrad = 0; % 1 for numerical gradient, 0 for analytical
end

if isfield(EstimOpt,'HessEstFix') == 0 || (EstimOpt.HessEstFix ~= 0 && EstimOpt.HessEstFix ~= 1 && EstimOpt.HessEstFix ~= 2 && EstimOpt.HessEstFix ~= 3&& EstimOpt.HessEstFix ~= 4)
    EstimOpt.HessEstFix = 0; % 0 = use optimization Hessian, 1 = use jacobian-based (BHHH) Hessian, 2 - use high-precision jacobian-based (BHHH) Hessian 3 - use numerical Hessian
end

if isfield(EstimOpt,'ApproxHess') == 0 || (EstimOpt.ApproxHess ~= 0 && EstimOpt.ApproxHess ~= 1)
    EstimOpt.ApproxHess = 1;
end

if isfield(EstimOpt,'NBtype') == 0 || isempty(EstimOpt.NBtype)
    EstimOpt.NBtype = 2; % 2 = NB2 (default), 1 = NB1
end

if isfield(EstimOpt,'RealMin') == 0 || (EstimOpt.RealMin ~= 0 && EstimOpt.RealMin ~= 1)
    EstimOpt.RealMin = 0;
end

%% OptimOpt

if isfield(EstimOpt, 'ConstVarActive') == 0 || EstimOpt.ConstVarActive == 0 % no contstaints on parameters
    OptimOpt = optimoptions('fminunc');
    OptimOpt.Algorithm = 'quasi-newton'; %'trust-region';
elseif EstimOpt.ConstVarActive == 1 % there are some constraints on parameters
    OptimOpt = optimoptions('fmincon');
    OptimOpt.Algorithm = 'interior-point'; %'sqp'; 'active-set'; 'trust-region-reflective';
end


OptimOpt.GradObj = 'on'; %'off';
% OptimOpt.FinDiffType = 'central'; % ('forward')
% OptimOpt.Hessian = 'user-supplied'; % ('off'), only used by trust-region

if isfield(EstimOpt,'FunctionTolerance')
    OptimOpt.FunctionTolerance = EstimOpt.FunctionTolerance; % df / gradient precision level
elseif isfield(EstimOpt,'eps')
    OptimOpt.FunctionTolerance = EstimOpt.eps;
end
if isfield(EstimOpt,'StepTolerance')
    OptimOpt.StepTolerance = EstimOpt.TolX; % step precision level
elseif isfield(EstimOpt,'eps')
    OptimOpt.StepTolerance = EstimOpt.eps;
end
if isfield(EstimOpt,'OptimalityTolerance')
    OptimOpt.OptimalityTolerance = EstimOpt.OptimalityTolerance; % dB precision level
elseif isfield(EstimOpt,'eps')
    OptimOpt.OptimalityTolerance = EstimOpt.eps;
end

OptimOpt.MaxIter = 1e4;
OptimOpt.FunValCheck = 'on';
OptimOpt.Diagnostics = 'off';
OptimOpt.MaxFunEvals = 1e5*size(INPUT.Xa,2); %Maximum number of function evaluations allowed (1000)
OptimOpt.OutputFcn = @outputf;


%% Estimate the constants-only model (same family as the model being estimated):

INPUT_0.Y = INPUT.Y;
INPUT_0.Xa = ones(size(INPUT.Y,1),1); % NCT*NP no longer holds once incomplete observations are removed
INPUT_0.MissingInd = INPUT.MissingInd;
INPUT_0.W = INPUT.W; %ones(EstimOpt.NP,1);
EstimOpt_0 = EstimOpt;
EstimOpt_0.ConstVarActive = 0;
EstimOpt_0.BActive = [];
EstimOpt_0.NVarA = 1;
EstimOpt_0.OPTIM = 1;
EstimOpt_0.Display = 0;
OptimOpt_0 = optimoptions('fminunc');
OptimOpt_0.Algorithm = 'quasi-newton';
OptimOpt_0.GradObj = 'off';
OptimOpt_0.Hessian = 'off';
OptimOpt_0.Display = 'off';
OptimOpt_0.FunValCheck= 'off';
OptimOpt_0.Diagnostics = 'off';
Results.CDM0 = CDM(INPUT_0,[],EstimOpt_0,OptimOpt_0);
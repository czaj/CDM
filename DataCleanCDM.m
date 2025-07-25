function [INPUT, Results, EstimOpt, OptimOpt] = DataCleanCDM(INPUT,EstimOpt)

% global TolB
% save tmp1
% return

inputnames = fieldnames(INPUT);
for i=1:length(inputnames)
    INPUT.(inputnames{(i)}) = double(INPUT.(inputnames{(i)}));
end

if isfield(INPUT,'MissingInd') == 0 || isempty(INPUT.MissingInd)
    INPUT.MissingInd = zeros(size(INPUT.Y));
end

EstimOpt.MissingCT = [];

if sum(INPUT.MissingInd) ~= 0
    MissingInd_tmp = reshape(INPUT.MissingInd,EstimOpt.NCT,EstimOpt.NP);
    MissingCT = MissingInd_tmp == 1; % missing NCT
    MissingP = sum(MissingCT,1) == EstimOpt.NCT; % respondents with all NCT missing
    
    if sum(MissingP) > 0 % respondents with 0 NCTs - remove from INPUT
        MissingPrep = reshape(MissingP(ones(EstimOpt.NCT,1),:),EstimOpt.NCT*EstimOpt.NP,1);
        INPUT_fields = fields(INPUT);
        for i = 1:size(INPUT_fields,1)
            tmp = INPUT.(INPUT_fields{i});
                tmp(MissingPrep,:) = [];
                INPUT.(INPUT_fields{i}) = tmp;
        end
        cprintf(rgb('DarkOrange'), ['WARNING: Dataset includes ', num2str(sum(MissingP)), ' respondents with 0 completed choice tasks. Adjusting NP from ', num2str(EstimOpt.NP), ' to ',num2str(EstimOpt.NP-sum(MissingP)) ,'.\n'])
        EstimOpt.NP = EstimOpt.NP - sum(MissingP);
    end
    
    
    if any(MissingCT(:)) > 0 % respondents with missing NCT - replace Xa and Y with NaN
        %cprintf ('text', 'The dataset contains %d choice tasks with missing responses (out of the total of %d choice tasks).\n', sum(sum(MissingCT)),numel(MissingCT))
        cprintf ('text', ['The dataset contains ',num2str(sum(sum(MissingCT))),' choice tasks with missing responses (out of the total of ',num2str(numel(MissingCT)) ,' choice tasks).\n'])
        INPUT.Y(INPUT.MissingInd == 1) = NaN;
        INPUT.Xa(INPUT.MissingInd == 1,:) = NaN;

    end
    EstimOpt.NCTMiss = EstimOpt.NCT - sum(MissingCT,1)';
else
    EstimOpt.NCTMiss = EstimOpt.NCT*ones(EstimOpt.NP,1);
end

EstimOpt.NObs = sum(EstimOpt.NCTMiss);

if isfield(INPUT,'W') && ~isempty(INPUT.W)
    if size(INPUT.W(:)) ~= size(INPUT.Y(:))
        error('Incorrect size of the weights vector')
    else
        INPUT.W = INPUT.W(:);
        INPUT.W = INPUT.W(1:EstimOpt.NCT:end);
        if sum(INPUT.W) ~= EstimOpt.NP
            cprintf(rgb('DarkOrange'), ['WARNING: Scaling weights for unit mean. \n'])
            INPUT.W = INPUT.W - mean(INPUT.W) + 1;
        end
    end
else
    INPUT.W = ones(EstimOpt.NP,1);
end

if isfield(EstimOpt,'RobustStd') == 0
    EstimOpt.RobustStd = 0; % do not use robust standard errors
end

EstimOpt.NVarA = size(INPUT.Xa,2); % Number of attributes

%EstimOpt.NVarZinf = size(INPUT.Xzinf,2);

EstimOpt.Display = 1;
if isfield(EstimOpt,'HaltonSkip') == 0
    EstimOpt.HaltonSkip = 1; % specify no of rows in halton sequence to skip (default=1)
end
if isfield(EstimOpt,'HaltonLeap') == 0
    EstimOpt.HaltonLeap = 0; % specify no of rows in halton sequence to leap (default=0)
end

if isfield(EstimOpt,'Draws') == 0
    EstimOpt.Draws = 6; % specify draws type (default = Sobol with scrambling)
end

if isfield(EstimOpt,'NRep') == 0
    EstimOpt.NRep = 1e3; % specify no. of draws
end

EstimOpt.Seed1 = 179424673;
EstimOpt.Seed2 = 7521436817;

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

if isfield(EstimOpt,'RealMin') == 0 || (EstimOpt.RealMin ~= 0 && EstimOpt.RealMin ~= 1)
    EstimOpt.RealMin = 0;
end

EstimOpt.NSdSim = 1e4; % number of draws for simulating standard deviations


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


%% Estimate constants-only Poisson model:

INPUT_0.Y = INPUT.Y;
INPUT_0.Xa = ones(EstimOpt.NCT*EstimOpt.NP,1);
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
Results.POISS0 = CDM(INPUT_0,[],EstimOpt_0,OptimOpt_0);
function nFail = test_CDM_recovery(verbose)
% test_CDM_recovery estimates each model on simulated data with known parameters
% and checks that the true values fall inside a 4 s.e. band around the estimates.
%
% Returns the number of failures.

if nargin < 1 || isempty(verbose)
    verbose = 1;
end
cdm_test_paths;

global B_backup %#ok<GVMIS>

nFail = 0;
nPass = 0;
nSE = 4;      % width of the acceptance band, in standard errors
n = 5000;

if verbose
    fprintf('test_CDM_recovery (n = %d, band = %d s.e.)\n',n,nSE);
    fprintf('  %-26s %10s %10s %10s   %s\n','model','max |t|','LL','k','status');
end

rng(987654321)
Xa = [ones(n,1),unifrnd(-1,1,n,1),unifrnd(0,2,n,1)];
beta = [1;-0.7;0.4];
Lam = exp(Xa*beta);
Xzinf = [ones(n,1),unifrnd(-1,1,n,1)];
delta = [0.3;-0.2];
theta = 2;

cases = {};

% Poisson
cases{end+1} = struct('name','Poisson', ...
    'Y',poissrnd(Lam),'EstimOpt',struct('NB',0,'Zinf',0), ...
    'Xzinf',[],'btrue',beta);

% Negative binomial: theta = exp(-gamma) with a constant-only Xnb
cases{end+1} = struct('name','Negative Binomial', ...
    'Y',nbinrnd(theta,theta./(theta+Lam)),'EstimOpt',struct('NB',1,'Zinf',0), ...
    'Xzinf',[],'btrue',[beta;-log(theta)]);

% Zero-truncated Poisson: resample the zeros so the draw is conditional on Y > 0
Ytrunc = poissrnd(Lam);
while any(Ytrunc == 0)
    idx = Ytrunc == 0;
    Ytrunc(idx) = poissrnd(Lam(idx));
end
cases{end+1} = struct('name','Zero-truncated Poisson', ...
    'Y',Ytrunc,'EstimOpt',struct('NB',0,'Zinf',0,'Truncated',1), ...
    'Xzinf',[],'btrue',beta);

% Zero-inflated Poisson: (1 - q) is the probability of being in the count regime
q = 1./(1+exp(-Xzinf*delta));
cases{end+1} = struct('name','Zero-inflated Poisson', ...
    'Y',poissrnd(Lam).*binornd(1,1-q),'EstimOpt',struct('NB',0,'Zinf',1), ...
    'Xzinf',Xzinf,'btrue',[beta;delta]);

for c = 1:numel(cases)
    tc = cases{c};
    B_backup = [];
    status = 'ok';
    maxT = NaN; LL = NaN; k = NaN;

    try
        clear INPUT EstimOpt
        INPUT.Xa = Xa;
        INPUT.Y = tc.Y;
        INPUT.MissingInd = zeros(n,1);
        if ~isempty(tc.Xzinf)
            INPUT.Xzinf = tc.Xzinf;
        end

        EstimOpt = tc.EstimOpt;
        EstimOpt.NP = n;
        EstimOpt.NCT = 1;
        EstimOpt.Display = 0;
        EstimOpt.NamesA = {'Cons';'Var1';'Var2'};

        [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
        OO.Display = 'off';
        OO.OutputFcn = [];
        OO.Diagnostics = 'off';

        R = CDM(IN,Res,EO,OO);

        LL = R.LL;
        k = numel(R.bhat);

        if k ~= numel(tc.btrue)
            status = sprintf('FAIL: estimated %d parameters, expected %d',k,numel(tc.btrue));
        elseif ~all(isfinite(R.std))
            status = 'FAIL: non-finite standard errors';
        else
            maxT = max(abs(R.bhat - tc.btrue)./R.std);
            if maxT > nSE
                [~,w] = max(abs(R.bhat - tc.btrue)./R.std);
                status = sprintf('FAIL: parameter %d off by %.2f s.e.',w,maxT);
            end
        end
    catch ME
        status = ['FAIL: ' ME.message];
    end

    if strcmp(status,'ok')
        nPass = nPass + 1;
    else
        nFail = nFail + 1;
    end

    if verbose
        fprintf('  %-26s %10.2f %10.1f %10d   %s\n',tc.name,maxT,LL,k,status);
    end
end

if verbose
    fprintf('test_CDM_recovery: %d passed, %d failed\n',nPass,nFail);
    if nFail == 0
        disp('test_CDM_recovery passed')
    end
end
end

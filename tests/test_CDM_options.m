function nFail = test_CDM_options(verbose)
% test_CDM_options exercises the structural paths of CDM/DataCleanCDM that are
% not covered by the likelihood grid: parameter indexing, Hessian options,
% Display, overdispersion covariates, weights and missing data.
%
% Returns the number of unexpected failures. Checks tagged KNOWN document open
% issues from the audit and are reported without failing the suite.

if nargin < 1 || isempty(verbose)
    verbose = 1;
end
cdm_test_paths;

global B_backup %#ok<GVMIS>
B_backup = [];

nFail = 0;
nKnown = 0;
nPass = 0;

    function report(id,name,ok,known,detail)
        if ok
            nPass = nPass + 1; tag = 'ok';
            if known, tag = ['ok (' id ' appears fixed - update KNOWN)']; end
        elseif known
            nKnown = nKnown + 1; tag = ['KNOWN ' id];
        else
            nFail = nFail + 1; tag = ['FAIL ' id];
        end
        if verbose
            fprintf('  %-10s %-52s %s\n',tag,name,detail);
        end
    end

if verbose
    fprintf('test_CDM_options\n');
end

rng(31415)
n = 1200;
Xa = [ones(n,1),unifrnd(-1,1,n,1)];
beta = [1;-0.6];
Lam = exp(Xa*beta);


%% B3 - the count equation must read only b(1:NVarA)
% Trailing parameters belong to other blocks (overdispersion, zero inflation,
% and anything added later). Appending one must not move the likelihood.
ok = true; detail = '';
try
    for Truncated = [0 1 2]
        Y = poissrnd(Lam);
        if Truncated > 0, Y = max(Y,1); end
        Y = min(Y,6);
        E = struct('NB',0,'Zinf',0,'Truncated',Truncated,'Censored',6, ...
            'NVarA',2,'NVarNB',0,'NVarZinf',0,'RealMin',0);
        f1 = LL_CDM(Y,Xa,zeros(n,0),zeros(n,0),E,beta);
        f2 = LL_CDM(Y,Xa,zeros(n,0),zeros(n,0),E,[beta;99]);
        if ~isequal(f1,f2)
            ok = false;
            detail = sprintf('Truncated=%d: LL moved by %.4g',Truncated,max(abs(f1-f2)));
            break
        end
    end
catch ME
    ok = false; detail = ME.message;
end
report('B3','count equation ignores trailing parameters',ok,false,detail);


%% B1 - RealMin must clamp at realmin, not at the flag value
ok = true; detail = '';
try
    Y = min(max(nbinrnd(2,2./(2+Lam)),1),6);
    E = struct('NB',1,'Zinf',0,'Truncated',1,'Censored',6, ...
        'NVarA',2,'NVarNB',1,'NVarZinf',0,'RealMin',0);
    f0 = LL_CDM(Y,Xa,ones(n,1),zeros(n,0),E,[beta;0.2]);
    E.RealMin = 1;
    f1 = LL_CDM(Y,Xa,ones(n,1),zeros(n,0),E,[beta;0.2]);
    ok = max(abs(f0-f1)) < 1e-12;
    detail = sprintf('max |dLL| between RealMin 0/1 = %.3g',max(abs(f0-f1)));
catch ME
    ok = false; detail = ME.message;
end
report('B1','zero-truncated censored NB invariant to RealMin',ok,false,detail);


%% B10 - DataCleanCDM must not override EstimOpt.Display
ok = true; detail = '';
try
    clear INPUT EstimOpt
    INPUT.Xa = Xa; INPUT.Y = poissrnd(Lam); INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0;
    [~,~,EO] = DataCleanCDM(INPUT,EstimOpt);
    ok = EO.Display == 0;
    detail = sprintf('Display after DataCleanCDM = %d',EO.Display);
catch ME
    ok = false; detail = ME.message;
end
report('B10','EstimOpt.Display = 0 is preserved',ok,false,detail);


%% B6 - HessEstFix = 4 must not crash
ok = true; detail = '';
try
    clear INPUT EstimOpt
    INPUT.Xa = Xa; INPUT.Y = poissrnd(Lam); INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0; EstimOpt.HessEstFix = 4;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO = quietOptim(OO);
    R = CDM(IN,Res,EO,OO);
    ok = all(isfinite(R.std));
    detail = sprintf('max s.e. = %.4f',max(R.std));
catch ME
    ok = false; detail = ME.message;
end
report('B6','HessEstFix = 4 falls back to a numerical Hessian',ok,false,detail);


%% B4 - covariates of overdispersion
ok = true; detail = '';
try
    clear INPUT EstimOpt
    INPUT.Xa = Xa; INPUT.Y = nbinrnd(2,2./(2+Lam));
    INPUT.Xnb = unifrnd(-1,1,n,1);
    INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.NB = 1; EstimOpt.Display = 0;
    EstimOpt.NamesA = {'Cons';'Var1'};
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO = quietOptim(OO);
    R = CDM(IN,Res,EO,OO);
    % NVarNB is resolved inside CDM (it prepends a constant to INPUT.Xnb)
    nExpected = R.EstimOpt.NVarA + R.EstimOpt.NVarNB;
    ok = numel(R.bhat) == nExpected && numel(R.EstimOpt.NamesNB) == R.EstimOpt.NVarNB ...
        && all(isfinite(R.std));
    detail = sprintf('%d of %d parameters, %d labels',numel(R.bhat),nExpected,numel(R.EstimOpt.NamesNB));
catch ME
    ok = false; detail = ME.message;
end
report('B4','NVarNB > 1 (covariates of overdispersion)',ok,false,detail);


%% B7 - respondent weights with NCT > 1
ok = true; detail = '';
try
    clear INPUT EstimOpt
    nn = 300; nct = 4; N = nn*nct;
    INPUT.Xa = [ones(N,1),unifrnd(-1,1,N,1)];
    INPUT.Y = poissrnd(exp(INPUT.Xa*beta));
    INPUT.MissingInd = zeros(N,1);
    INPUT.W = repelem(unifrnd(0.5,1.5,nn,1),nct,1);
    EstimOpt.NP = nn; EstimOpt.NCT = nct; EstimOpt.Display = 0;
    base.Xa = INPUT.Xa; base.Y = INPUT.Y; base.MissingInd = INPUT.MissingInd;

    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    Rw = CDM(IN,Res,EO,quietOptim(OO));

    % no weights at all
    [IN,Res,EO,OO] = DataCleanCDM(base,EstimOpt);
    Rn = CDM(IN,Res,EO,quietOptim(OO));

    % constant weights must reproduce the unweighted fit exactly
    INPUTc = base; INPUTc.W = 3*ones(N,1);
    [IN,Res,EO,OO] = DataCleanCDM(INPUTc,EstimOpt);
    Rc = CDM(IN,Res,EO,quietOptim(OO));

    reachedLik = ~all(Rw.INPUT.W == 1);
    constSame = max(abs(Rc.bhat - Rn.bhat)) < 1e-6;
    weightsBite = max(abs(Rw.bhat - Rn.bhat)) > 1e-4;

    ok = reachedLik && constSame && weightsBite;
    detail = sprintf('reached likelihood %d, constant weights neutral %d (d=%.1g), varying weights move bhat %d (d=%.3g)', ...
        reachedLik,constSame,max(abs(Rc.bhat-Rn.bhat)),weightsBite,max(abs(Rw.bhat-Rn.bhat)));
catch ME
    ok = false; detail = ME.message;
end
report('B7','weights survive DataCleanCDM -> CDM with NCT > 1',ok,false,detail);


%% B5 - partially missing panel observations
ok = true; detail = '';
try
    clear INPUT EstimOpt
    nn = 300; nct = 4; N = nn*nct;
    INPUT.Xa = [ones(N,1),unifrnd(-1,1,N,1)];
    INPUT.Y = poissrnd(exp(INPUT.Xa*beta));
    INPUT.MissingInd = zeros(N,1); INPUT.MissingInd(2:7:N) = 1;
    EstimOpt.NP = nn; EstimOpt.NCT = nct; EstimOpt.Display = 0;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO = quietOptim(OO);
    R = CDM(IN,Res,EO,OO);
    ok = isfinite(R.LL) && ~any(isnan(IN.Y));
    detail = sprintf('LL = %.4f, NaN rows left = %d',R.LL,sum(isnan(IN.Y)));
catch ME
    ok = false; detail = ME.message;
end
report('B5','missing observations removed when NCT > 1',ok,false,detail);


%% B9 - the dependent variable must be validated, not silently mis-handled
ok = true; detail = '';
try
    clear INPUT EstimOpt
    Yraw = poissrnd(Lam);
    nAbove = sum(Yraw > 6);
    base.Xa = Xa; base.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0; EstimOpt.Censored = 6;

    % (a) Y above the censoring point is censored by CDM, not treated as uncensored
    INPUT = base; INPUT.Y = Yraw;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    Rraw = CDM(IN,Res,EO,quietOptim(OO));

    INPUT = base; INPUT.Y = min(Yraw,6);
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    Rcens = CDM(IN,Res,EO,quietOptim(OO));

    okA = abs(Rraw.LL - Rcens.LL) < 1e-6;

    % (b) non-integer and negative counts are rejected
    okB = errorsOut(@() runCDM(base,Yraw + 0.5,EstimOpt));
    okC = errorsOut(@() runCDM(base,Yraw - 1,EstimOpt));

    % (c) a zero-truncated model with zeros in the data is rejected
    EstimOptT = EstimOpt; EstimOptT.Truncated = 1; EstimOptT.Censored = 0;
    okD = errorsOut(@() runCDM(base,Yraw,EstimOptT));

    ok = okA && okB && okC && okD;
    detail = sprintf('%d obs above the cut censored (dLL = %.2g); rejects non-integer %d, negative %d, zeros under truncation %d', ...
        nAbove,Rraw.LL-Rcens.LL,okB,okC,okD);
catch ME
    ok = false; detail = ME.message;
end
report('B9','Y is validated and censored at EstimOpt.Censored',ok,false,detail);


%% B8 / B16 - the likelihood must stay finite in the corners of the parameter space
ok = true; detail = '';
try
    worst = 0;
    % (a) near-Poisson NB: theta = exp(-ln alpha) overflows gamma() past ln alpha = -5.2
    Yc = min(poissrnd(Lam),6);
    for lnAlpha = [-3 -6 -9]
        E = struct('NB',1,'Zinf',0,'Truncated',0,'Censored',6, ...
            'NVarA',2,'NVarNB',1,'NVarZinf',0,'RealMin',0);
        [f,g] = LL_CDM(Yc,Xa,ones(n,1),zeros(n,0),E,[beta;lnAlpha]);
        if ~all(isfinite(f)) || ~all(isfinite(g(:)))
            ok = false;
            detail = sprintf('non-finite LL/gradient at ln(alpha) = %d',lnAlpha);
            break
        end
        worst = max(worst,max(abs(g(:))));
    end

    % (b) zero-truncated Poisson at a tiny lambda: 1-exp(-lambda) cancels catastrophically
    if ok
        Xsmall = [ones(n,1),zeros(n,1)];
        E = struct('NB',0,'Zinf',0,'Truncated',1,'Censored',0, ...
            'NVarA',2,'NVarNB',0,'NVarZinf',0,'RealMin',0);
        for lnLam = [-8 -16 -30]
            [f,g] = LL_CDM(ones(n,1),Xsmall,zeros(n,0),zeros(n,0),E,[lnLam;0]);
            if ~all(isfinite(f)) || ~all(isfinite(g(:)))
                ok = false;
                detail = sprintf('non-finite LL/gradient at log(lambda) = %d',lnLam);
                break
            end
        end
    end

    if ok
        detail = sprintf('finite down to ln(alpha) = -9 and log(lambda) = -30 (max |g| = %.3g)',worst);
    end
catch ME
    ok = false; detail = ME.message;
end
report('B8/B16','likelihood finite at extreme theta and lambda',ok,false,detail);


if verbose
    fprintf('test_CDM_options: %d passed, %d known issues, %d failed\n',nPass,nKnown,nFail);
    if nFail == 0
        disp('test_CDM_options passed')
    end
end
end


function OptimOpt = quietOptim(OptimOpt)
OptimOpt.Display = 'off';
OptimOpt.OutputFcn = [];
OptimOpt.Diagnostics = 'off';
end


function runCDM(base,Y,EstimOpt)
% runCDM estimates one model - used to check that bad input is rejected.
INPUT = base;
INPUT.Y = Y;
[IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
CDM(IN,Res,EO,quietOptim(OO));
end


function tf = errorsOut(fun)
% errorsOut is true when fun() raises an error.
tf = false;
try
    fun();
catch
    tf = true;
end
end

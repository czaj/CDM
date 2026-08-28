function nFail = test_CDM_post(verbose)
% test_CDM_post checks the post-estimation results: the expected-trip formulas
% behind consumer surplus, the recovery of consumer surplus on simulated data,
% and the model comparison tests.
%
% Returns the number of failures.

if nargin < 1 || isempty(verbose)
    verbose = 1;
end
cdm_test_paths;

global B_backup %#ok<GVMIS>

nFail = 0;
nPass = 0;

    function report(name,ok,detail)
        if ok
            nPass = nPass + 1; tag = 'ok';
        else
            nFail = nFail + 1; tag = 'FAIL';
        end
        if verbose
            fprintf('  %-6s %-46s %s\n',tag,name,detail);
        end
    end

if verbose
    fprintf('test_CDM_post\n');
end

%% CDM_EY against Monte Carlo, one variant at a time
% Constant lambda (and theta) so the sample mean of a large draw is the check.
rng(20260829)
m = 4e5;
lam = 2.3; th = 1.7;
Xa1 = ones(1,1); Xnb1 = ones(1,1);
E = struct('NVarA',1,'NVarNB',0,'NVarZinf',0,'NB',0,'Zinf',0,'Truncated',0);

    function checkEY(name,EstimOptCase,bvec,draw)
        Ehat = CDM_EY(bvec,Xa1,Xnb1,zeros(1,0),EstimOptCase);
        mc = mean(draw);
        rel = abs(Ehat-mc)/mc;
        report(name,rel < 0.01,sprintf('formula %.4f vs Monte Carlo %.4f (%.2f%%)',Ehat,mc,100*rel));
    end

% Poisson, no truncation
checkEY('E[Y] Poisson',E,log(lam),poissrnd(lam,m,1));

% Poisson, zero-truncated
Ep = E; Ep.Truncated = 1;
yt = poissrnd(lam,m,1);
while any(yt == 0), yt(yt == 0) = poissrnd(lam,sum(yt == 0),1); end
checkEY('E[Y|Y>0] zero-truncated Poisson',Ep,log(lam),yt);

% Poisson, endogenous stratification: Y-1 is Poisson
Ep2 = E; Ep2.Truncated = 2;
checkEY('E*[Y] stratified Poisson',Ep2,log(lam),poissrnd(lam,m,1)+1);

% Negative Binomial, no truncation
En = E; En.NB = 1; En.NVarNB = 1;
bn = [log(lam);-log(th)]; % theta = exp(-gamma)
checkEY('E[Y] Negative Binomial',En,bn,nbinrnd(th,th/(th+lam),m,1));

% NB, zero-truncated
En1 = En; En1.Truncated = 1;
yn = nbinrnd(th,th/(th+lam),m,1);
while any(yn == 0), yn(yn == 0) = nbinrnd(th,th/(th+lam),sum(yn == 0),1); end
checkEY('E[Y|Y>0] zero-truncated NB',En1,bn,yn);

% NB, endogenous stratification: size-biased, drawn by resampling with weight y
En2 = En; En2.Truncated = 2;
yb = nbinrnd(th,th/(th+lam),m,1);
idx = randsample(m,m,true,yb);
checkEY('E*[Y] stratified NB',En2,bn,yb(idx));

%% Consumer surplus recovered on simulated data
rng(777)
n = 4000;
Xa = [ones(n,1),unifrnd(5,50,n,1),unifrnd(-1,1,n,1)];
beta = [1.8;-0.05;0.3];
Lam = exp(Xa*beta);
csTrue = -1/beta(2);

ok = true; detail = '';
try
    B_backup = [];
    clear INPUT EstimOpt
    INPUT.Xa = Xa; INPUT.Y = poissrnd(Lam); INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0;
    EstimOpt.NamesA = {'Cons';'Travel cost';'Income'};
    EstimOpt.TC_par_index = 2;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO.Display = 'off'; OO.OutputFcn = [];
    R = CDM(IN,Res,EO,OO);
    R.EstimOpt.Display = 0;
    R = CDM_post(R,struct('NSim',300,'Display',0));

    eyTrue = mean(Lam);
    inTrip = csTrue >= R.CS.PerTrip(3) && csTrue <= R.CS.PerTrip(5);
    inPers = eyTrue*csTrue >= R.CS.PerPerson(3) && eyTrue*csTrue <= R.CS.PerPerson(5);
    inEY = eyTrue >= R.CS.EYsim(1) && eyTrue <= R.CS.EYsim(3);
    ok = inTrip && inPers && inEY;
    detail = sprintf('CS/trip %.2f (true %.2f), CS/person %.2f (true %.2f), E[Y] %.3f (true %.3f)', ...
        R.CS.PerTrip(1),csTrue,R.CS.PerPerson(1),eyTrue*csTrue,R.CS.EYmean,eyTrue);
    Rpoiss = R;
catch ME
    ok = false; detail = ME.message; Rpoiss = [];
end
report('consumer surplus inside the 95% interval',ok,detail);

% the predicted count distribution must integrate to one
ok = false; detail = '';
if ~isempty(Rpoiss)
    tot = sum(Rpoiss.CountFit.Predicted);
    ok = abs(tot-1) < 0.02;
    detail = sprintf('predicted frequencies sum to %.4f',tot);
end
report('predicted count distribution sums to one',ok,detail);

%% Model comparison
ok = true; detail = '';
try
    B_backup = [];
    clear INPUT EstimOpt
    th2 = 1.5;
    INPUT.Xa = Xa; INPUT.Y = nbinrnd(th2,th2./(th2+Lam)); INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0;
    EstimOpt.NamesA = {'Cons';'Travel cost';'Income'};

    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO.Display = 'off'; OO.OutputFcn = [];
    Rp = CDM(IN,Res,EO,OO);

    B_backup = [];
    EstimOpt.NB = 1;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO.Display = 'off'; OO.OutputFcn = [];
    Rn = CDM(IN,Res,EO,OO);

    T = CDM_compare(Rp,Rn,0);
    ok = T.pValue < 0.01 && contains(T.type,'boundary-corrected');
    detail = sprintf('LR %.1f, p %.2g on overdispersed data',T.statistic,T.pValue);
catch ME
    ok = false; detail = ME.message;
end
report('LR test detects overdispersion',ok,detail);

% Voung: on data with a genuine excess of zeros, the zero-inflated model must win
ok = true; detail = '';
try
    B_backup = [];
    clear INPUT EstimOpt
    Xz = [ones(n,1),unifrnd(-1,1,n,1)];
    q = 1./(1+exp(-Xz*[0.3;-0.8]));
    INPUT.Xa = Xa; INPUT.Y = poissrnd(Lam).*binornd(1,1-q); INPUT.MissingInd = zeros(n,1);
    EstimOpt.NP = n; EstimOpt.NCT = 1; EstimOpt.Display = 0;
    EstimOpt.NamesA = {'Cons';'Travel cost';'Income'};

    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO.Display = 'off'; OO.OutputFcn = [];
    Rp = CDM(IN,Res,EO,OO);

    B_backup = [];
    INPUT.Xzinf = Xz; EstimOpt.Zinf = 1;
    [IN,Res,EO,OO] = DataCleanCDM(INPUT,EstimOpt);
    OO.Display = 'off'; OO.OutputFcn = [];
    Rz = CDM(IN,Res,EO,OO);

    evalc('T = CDM_compare(Rp,Rz,0);'); % voung_test prints unconditionally
    ok = isscalar(T.statistic) && T.statistic < -1.96 && contains(T.conclusion,'zero-inflated');
    detail = sprintf('Voung %.2f -> %s',T.statistic,T.conclusion);
catch ME
    ok = false; detail = ME.message;
end
report('Voung test detects zero inflation',ok,detail);

if verbose
    fprintf('test_CDM_post: %d passed, %d failed\n',nPass,nFail);
    if nFail == 0
        disp('test_CDM_post passed')
    end
end
end

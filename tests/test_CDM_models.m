function nFail = test_CDM_models(verbose)
% test_CDM_models checks the model variants added on top of the Poisson/NB core:
% NB1, the hurdle model, and the travel cost helper.
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
            fprintf('  %-6s %-44s %s\n',tag,name,detail);
        end
    end

if verbose
    fprintf('test_CDM_models\n');
end

rng(20260829)
n = 5000;
Xa = [ones(n,1),unifrnd(5,50,n,1),unifrnd(-1,1,n,1)];
beta = [1.8;-0.05;0.3];
Lam = exp(Xa*beta);

    function [R,err] = fitCDM(Y,extraOpt,Xzinf)
        R = []; err = '';
        try
            B_backup = [];
            IN = struct('Xa',Xa,'Y',Y,'MissingInd',zeros(n,1));
            if nargin > 2 && ~isempty(Xzinf), IN.Xzinf = Xzinf; end
            EO = struct('NP',n,'NCT',1,'Display',0, ...
                'NamesA',{{'Cons';'Travel cost';'Income'}},'TC_par_index',2);
            f = fieldnames(extraOpt);
            for q = 1:numel(f), EO.(f{q}) = extraOpt.(f{q}); end
            [IN,Res,EO,OO] = DataCleanCDM(IN,EO);
            OO.Display = 'off'; OO.OutputFcn = [];
            R = CDM(IN,Res,EO,OO);
        catch ME
            err = ME.message;
        end
    end

%% NB1 recovers its parameters and is distinguishable from NB2
alphaTrue = 0.8;
thNB1 = Lam/alphaTrue; % NB1: variance lambda*(1+alpha)
Y1 = nbinrnd(thNB1,thNB1./(thNB1+Lam));

[R1,err] = fitCDM(Y1,struct('NB',1,'NBtype',1));
ok = isempty(err);
detail = err;
if ok
    tb = abs(R1.bhat(1:3)-beta)./R1.std(1:3);
    ta = abs(R1.bhat(4)-log(alphaTrue))/R1.std(4);
    ok = max([tb;ta]) < 4;
    detail = sprintf('beta %s, ln(alpha) %.3f (true %.3f), max |t| %.2f', ...
        mat2str(R1.bhat(1:3)',3),R1.bhat(4),log(alphaTrue),max([tb;ta]));
end
report('NB1 recovers beta and ln(alpha)',ok,detail);

% NB1 must beat NB2 on NB1 data
[R2,err2] = fitCDM(Y1,struct('NB',1,'NBtype',2));
ok = isempty(err) && isempty(err2) && R1.LL > R2.LL;
report('NB1 fits NB1 data better than NB2',ok, ...
    sprintf('LL NB1 %.1f vs NB2 %.1f',R1.LL,R2.LL));

% and the reverse on NB2 data
thNB2 = 2*ones(n,1);
Y2 = nbinrnd(thNB2,thNB2./(thNB2+Lam));
[R1b,~] = fitCDM(Y2,struct('NB',1,'NBtype',1));
[R2b,~] = fitCDM(Y2,struct('NB',1,'NBtype',2));
ok = ~isempty(R1b) && ~isempty(R2b) && R2b.LL > R1b.LL;
report('NB2 fits NB2 data better than NB1',ok, ...
    sprintf('LL NB2 %.1f vs NB1 %.1f',R2b.LL,R1b.LL));

%% NB-P estimates the power and nests NB1 and NB2
Ptrue = 1.5;
thP = Lam.^(2-Ptrue)/alphaTrue;
Yp = nbinrnd(thP,thP./(thP+Lam));

[Rp1,errp] = fitCDM(Yp,struct('NB',1,'NBtype',3));
ok = isempty(errp);
detail = errp;
if ok
    Phat = Rp1.bhat(end);
    ok = abs(Phat-Ptrue)/Rp1.std(end) < 4 && all(isfinite(Rp1.std));
    detail = sprintf('P %.3f (true %.3f), s.e. %.3f',Phat,Ptrue,Rp1.std(end));
end
report('NB-P recovers the power P',ok,detail);

% NB-P nests both, so on NB-P data it must fit at least as well as either
[Rp2,~] = fitCDM(Yp,struct('NB',1,'NBtype',2));
[Rp3,~] = fitCDM(Yp,struct('NB',1,'NBtype',1));
ok = ~isempty(Rp1) && ~isempty(Rp2) && ~isempty(Rp3) && ...
     Rp1.LL >= Rp2.LL - 1e-6 && Rp1.LL >= Rp3.LL - 1e-6;
report('NB-P fits at least as well as NB1 and NB2',ok, ...
    sprintf('LL NB-P %.1f, NB2 %.1f, NB1 %.1f',Rp1.LL,Rp2.LL,Rp3.LL));

% and on NB2 data the estimated power should sit near 2
[Rp4,~] = fitCDM(Y2,struct('NB',1,'NBtype',3));
ok = ~isempty(Rp4) && abs(Rp4.bhat(end)-2)/Rp4.std(end) < 4;
report('NB-P lands near P = 2 on NB2 data',ok, ...
    sprintf('P %.3f (s.e. %.3f)',Rp4.bhat(end),Rp4.std(end)));

% NB1 with censoring now has an analytical gradient, so it must not need NumGrad
[Rc,errc] = fitCDM(min(max(Y1,1),8),struct('NB',1,'NBtype',1,'Truncated',1,'Censored',8));
if isempty(Rc)
    ok = false; detail = errc;
else
    ok = Rc.EstimOpt.NumGrad == 0 && all(isfinite(Rc.std));
    detail = sprintf('NumGrad %d, max s.e. %.4f',Rc.EstimOpt.NumGrad,max(Rc.std));
end
report('NB1 with censoring keeps the analytical gradient',ok,detail);

%% Hurdle model
deltaTrue = [0.6;-0.03;0.5];
Ppart = 1./(1+exp(-Xa*deltaTrue));
part = binornd(1,Ppart);
Yt = poissrnd(Lam);
while any(Yt == 0), Yt(Yt == 0) = poissrnd(Lam(Yt == 0)); end
Yh = part.*Yt;

ok = true; detail = '';
Rh = [];
try
    B_backup = [];
    IN = struct('Xa',Xa,'Y',Yh,'MissingInd',zeros(n,1));
    EO = struct('NP',n,'NCT',1,'Display',0, ...
        'NamesA',{{'Cons';'Travel cost';'Income'}}, ...
        'NamesH',{{'Travel cost';'Income'}},'TC_par_index',2);
    [IN,Res,EO,OO] = DataCleanCDM(IN,EO);
    OO.Display = 'off'; OO.OutputFcn = [];
    EO.Display = 0;
    Rh = CDM_hurdle(IN,Res,EO,OO);

    tCount = abs(Rh.Count.bhat(1:3)-beta)./Rh.Count.std(1:3);
    tBin = abs(Rh.Binary.bhat-deltaTrue)./Rh.Binary.std;
    ok = max([tCount;tBin]) < 4;
    detail = sprintf('count max |t| %.2f, participation max |t| %.2f',max(tCount),max(tBin));
catch ME
    ok = false; detail = ME.message;
end
report('hurdle recovers both equations',ok,detail);

% the two parts must add up
ok = false; detail = '';
if ~isempty(Rh)
    ok = abs(Rh.LL-(Rh.Binary.LL+Rh.Count.LL)) < 1e-8 && ...
         abs(sum(Rh.LLdetailed)-Rh.LL) < 1e-6;
    detail = sprintf('LL %.4f = %.4f + %.4f, LLdetailed sums to %.4f', ...
        Rh.LL,Rh.Binary.LL,Rh.Count.LL,sum(Rh.LLdetailed));
end
report('hurdle log-likelihood separates cleanly',ok,detail);

% post-estimation on a hurdle: the zero cell must be predicted, not skipped
ok = false; detail = '';
if ~isempty(Rh)
    try
        Rh.EstimOpt.Display = 0;
        Rh = CDM_post(Rh,struct('NSim',200,'Display',0));
        obs0 = Rh.CountFit.Observed(1);
        pre0 = Rh.CountFit.Predicted(1);
        csTrue = -1/beta(2);
        ok = Rh.CountFit.k(1) == 0 && abs(obs0-pre0) < 0.02 && ...
             csTrue >= Rh.CS.PerTrip(3) && csTrue <= Rh.CS.PerTrip(5);
        detail = sprintf('P(Y=0) observed %.4f vs predicted %.4f, CS/trip %.2f (true %.2f)', ...
            obs0,pre0,Rh.CS.PerTrip(1),csTrue);
    catch ME
        detail = ME.message;
    end
end
report('hurdle post-estimation predicts the zeros',ok,detail);

%% Travel cost helper
ok = true; detail = '';
try
    Trip = struct();
    Trip.Distance = unifrnd(5,200,n,1);
    Trip.Time = Trip.Distance/70;
    Trip.CostPerKm = 0.35;
    Trip.Wage = unifrnd(15,60,n,1);
    Trip.ValuesTime = binornd(1,0.6,n,1);

    tcNone = TravelCost(Trip,struct('TimeRule','none','Display',0));
    tcFull = TravelCost(Trip,struct('TimeRule','full','Display',0));
    tcMix  = TravelCost(Trip,struct('TimeRule','full_or_zero','Display',0));
    tcFrac = TravelCost(Trip,struct('TimeRule','fraction','TimeFraction',[0 0.5 1],'Display',0));

    v = Trip.ValuesTime == 1;
    ok = all(tcNone <= tcMix & tcMix <= tcFull) ...
        && all(abs(tcMix(~v)-tcNone(~v)) < 1e-12) ...
        && all(abs(tcMix(v)-tcFull(v)) < 1e-12) ...
        && size(tcFrac,2) == 3 ...
        && all(abs(tcFrac(:,1)-tcNone) < 1e-12) ...
        && all(abs(tcFrac(:,3)-tcFull) < 1e-12) ...
        && abs(mean(Trip.Distance*2*0.35)-mean(tcNone)) < 1e-8;
    detail = sprintf('mean money %.2f, full wage %.2f, full-or-zero %.2f', ...
        mean(tcNone),mean(tcFull),mean(tcMix));
catch ME
    ok = false; detail = ME.message;
end
report('travel cost rules are internally consistent',ok,detail);

% the rule must be refused when the indicator it needs is missing
ok = false; detail = '';
try
    Trip2 = rmfield(Trip,'ValuesTime');
    try
        TravelCost(Trip2,struct('TimeRule','full_or_zero','Display',0));
        detail = 'accepted a full_or_zero call without Trip.ValuesTime';
    catch
        ok = true; detail = 'missing Trip.ValuesTime is rejected';
    end
catch ME
    detail = ME.message;
end
report('full_or_zero requires its indicator',ok,detail);

if verbose
    fprintf('test_CDM_models: %d passed, %d failed\n',nPass,nFail);
    if nFail == 0
        disp('test_CDM_models passed')
    end
end
end

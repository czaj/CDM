function nFail = test_CDM_hybrid_api(verbose)
% Standard preparation, both calling conventions and optimizer objective.
if nargin<1, verbose=1; end
cdm_test_paths;
nFail = 0;
checks = {@preparation,@callingOrders,@objectiveWrapper,@mappedStarts,@quietFit};
for j=1:numel(checks)
    try
        detail=checks{j}(); tag='ok';
    catch ME
        nFail=nFail+1; detail=ME.message; tag='FAIL';
    end
    if verbose, fprintf('  %-5s %-28s %s\n',tag,func2str(checks{j}),detail); end
end
if verbose, fprintf('test_CDM_hybrid_api: %d passed, %d failed\n',numel(checks)-nFail,nFail); end
end

function detail=preparation
[I,O]=example;
[I,R,O,P]=DataCleanCDM_hybrid(I,O);
assert(isempty(fieldnames(R)),'Preparation must not pretend to have optimized a null model.');
assert(numel(O.InitialValues)==O.Params && numel(O.ParameterNames)==O.Params);
assert(O.NP==numel(I.Y) && O.NVarMea==3 && O.NumGrad==0);
assert(strcmp(P.GradObj,'on') && O.SaveXLS==0 && O.SaveTxtOutput==0);
assert(isequal(I.W,(1:numel(I.Y))'/numel(I.Y)+0.5),'Weights were changed by preparation.');
assert(isnan(I.Xmea(3,2)),'Missing indicator was removed.');
[I2,~,O2]=DataCleanCDM_hybrid(I,O);
assert(isequaln(I,I2) && isequaln(O.InitialValues,O2.InitialValues),'Preparation is not idempotent.');
detail='four outputs, names/starts, unchanged rows/weights and no fictional null LL';
end

function detail=callingOrders
[I,O]=example;
R=CDM_hybrid(I,O);
P=optimset('Display','off','MaxIter',3);
Q=optimoptions('fminunc','Display','off','MaxIterations',3);
for optimizer={P,Q,[]}
    for prior={[],struct(),R}
        A=CDM_hybrid(I,prior{1},O,optimizer{1});
        B=CDM_hybrid(I,O,optimizer{1},prior{1});
        assert(isequal(A.bhat,B.bhat) && A.LL==B.LL,'Four-argument calling orders differ.');
        assert(isequal(A.bhat,R.bhat) && A.LL==R.LL);
    end
end
% Empty structs are valid optimizer options and model option containers.
A=CDM_hybrid(I,R,O,struct()); B=CDM_hybrid(I,O,struct(),R);
assert(A.LL==B.LL);
detail='standard/legacy calls agree with empty, struct and optimoptions optimizers';
end

function detail=objectiveWrapper
[I,O]=example; [I,~,O,P]=DataCleanCDM_hybrid(I,O);
b=O.InitialValues; b(O.Layout.delta)=0.18; b(O.Layout.gamma)=-0.1;
[lli,score]=LL_CDM_hybrid(b,I,O);
[f,g,H]=LL_CDM_hybrid_MATlike(I,O,P,b);
assert(f==-sum(lli) && isequal(g,-sum(score,1)') && isequal(H,score'*score));
numeric=zeros(size(b));
for j=1:numel(b)
    h=1e-5*(1+abs(b(j))); bp=b; bm=b; bp(j)=bp(j)+h; bm(j)=bm(j)-h;
    numeric(j)=(LL_CDM_hybrid_MATlike(I,O,P,bp)-LL_CDM_hybrid_MATlike(I,O,P,bm))/(2*h);
end
err=max(abs(g-numeric)./(1+abs(numeric)));
assert(err<1e-6,'Negative objective gradient mismatch %.3g.',err);
detail=sprintf('weighted sign/score/BHHH agreement; finite-difference error %.2g',err);
end

function detail=mappedStarts
[I,O]=example; O.LatentActive=false;
restricted=CDM_hybrid(I,O);
O.LatentActive=true;
R=CDM_hybrid(I,restricted,O,[]);
assert(R.bhat(R.Layout.delta)==0 && abs(R.LL-restricted.LL)<1e-10);
R2=CDM_hybrid(I,struct('Restricted',restricted),O,[]);
assert(abs(R2.LL-restricted.LL)<1e-10);
other=restricted; other.EstimOpt.NB=0;
try
    CDM_hybrid(I,other,O,[]);
    error('test:MissingRejection','Different count family accepted as a raw warm start.');
catch ME
    assert(strcmp(ME.identifier,'CDM_hybrid:Start'),'Unexpected family-mismatch error: %s',ME.message);
end
O.Start=R.bhat;
explicit=CDM_hybrid(I,other,O,[]);
assert(isequal(explicit.bhat,O.Start),'Explicit Start must override previous results.');
detail='restricted-to-joint names map; accumulator supported; count-family mismatch rejected';
end

function detail=quietFit
[I,O]=example; O.EvaluateOnly=false; O.ComputeHessian=false;
[I,Results,O,P]=DataCleanCDM_hybrid(I,O);
P.MaxIter=4; P.MaxFunEvals=100;
A=CDM_hybrid(I,Results,O,P); B=CDM_hybrid(I,O,P);
assert(abs(A.LL-B.LL)<1e-10 && max(abs(A.bhat-B.bhat))<1e-10);
assert(~isfield(A,'OutputXLSX') && strcmp(A.OptimOpt.GradObj,'on'));
detail='both entry forms follow identical analytic-gradient optimization without output';
end

function [I,O]=example
n=32; z=linspace(-1.5,1.5,n)';
I=struct('Y',1+mod((1:n)',5),'Xa',[ones(n,1),z], ...
    'Xstr',z,'Xmea',repmat(dec2bin(0:7)-'0'+1,4,1), ...
    'W',(1:n)'/n+0.5);
I.Xmea(3,2)=NaN;
O=struct('Display',0,'EvaluateOnly',true,'Truncated',1,'NB',1, ...
    'NNodes',20,'NamesA',{{'Constant';'Cost'}},'NamesStr',{{'Income'}}, ...
    'NamesMea',{{'Risk';'Concern';'Attention'}});
end

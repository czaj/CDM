function [nFail,outputDir] = test_CDM_hybrid_output(verbose)
% Standard Excel output, delta-method transformations and inference gates.
% Uses known covariance fixtures; no fit is run or evidence manufactured.
if nargin<1, verbose=1; end
cdm_test_paths;
outputDir = fullfile(tempdir,'CDM_hybrid_output_tests');
nFail = 0;
if ~ispc
    if verbose, fprintf('  SKIP standard XLSX output requires Windows Excel COM.\n'); end
    return
end
checks = {@naturalOutput,@inferenceGate,@standardCompatibility};
for k=1:numel(checks)
    try
        detail=checks{k}(outputDir); tag='ok';
    catch ME
        nFail=nFail+1; detail=getReport(ME,'basic','hyperlinks','off'); tag='FAIL';
    end
    if verbose, fprintf('  %-5s %-25s %s\n',tag,func2str(checks{k}),detail); end
end
if verbose, fprintf('test_CDM_hybrid_output: %d passed, %d failed\n',numel(checks)-nFail,nFail); end
end

function detail=naturalOutput(outDir)
[R,O]=fixture(outDir);
V=R.ihess; b=R.bhat; L=R.Layout;
R=CDM_hybrid_output(R,O);
assert(R.OutputInferenceValid && exist(R.OutputXLSX,'file')==2,'Missing validated output.');
assert(max(abs(R.DetailsRaw(:,1)-b))==0,'Raw parameter vector changed.');
a=exp(b(L.logalpha));
assert(abs(R.DetailsNB(1,1)-a)<1e-12 && abs(R.DetailsNB(1,3)-a*sqrt(V(L.logalpha,L.logalpha)))<1e-12);
assert(isnan(R.DetailsNB(1,4)),'Alpha has an inappropriate zero-test p-value.');
ix=L.loading(1);
assert(abs(R.DetailsMea1(1,1)-exp(b(ix)))<1e-12);
assert(abs(R.DetailsMea1(1,3)-exp(b(ix))*sqrt(V(ix,ix)))<1e-12);
assert(isnan(R.DetailsMea1(1,4)),'Anchored positive loading has a zero-test p-value.');
ix=L.thresholds{1}; gap=exp(b(ix(2)));
expectedVariance=V(ix(1),ix(1))+gap^2*V(ix(2),ix(2))+2*gap*V(ix(1),ix(2));
assert(abs(R.DetailsMea1(3,1)-(b(ix(1))+gap))<1e-12);
assert(abs(R.DetailsMea1(3,3)-sqrt(expectedVariance))<1e-12,'Cutoff SE omits covariance.');
assert(numel(R.OutputBlocks)==8,'Missing count, dispersion, LV, structural, item or raw block.');
assert(isequal(cellstr(sheetnames(R.OutputXLSX)),{'Results';'Script';'Log';'Word 2dp';'Word 4dp'}),'Wrong standard workbook tabs.');
C=readcell(R.OutputXLSX,'Sheet','Results');
assert(any(strcmp(C(:,1),'Travel cost unit')),'Application scaling note missing.');
row=find(strcmp(C(:,1),'alpha'),1);
assert(abs(C{row,3}-a)<1e-12 && abs(C{row,5}-a*sqrt(V(L.logalpha,L.logalpha)))<1e-12,'Excel numeric values differ.');
assert(isnumeric(C{row,3}) && isnumeric(C{row,5}),'XLSX estimates must remain numeric.');
W=readcell(R.OutputXLSX,'Sheet','Word 4dp');
assert(size(W,1)>numel(b),'Word-friendly output is incomplete.');
detail='natural alpha/loadings/cutoffs, covariance terms, raw vector and 5 standard tabs';
end

function detail=inferenceGate(outDir)
[R,O]=fixture(outDir); O.OutputFilePrefix='hybrid_unchecked';
R.QuadratureStable=false;
% Export controls cannot override statistical validity stored in the fit.
O.InferenceValid=true; O.QuadratureStable=true;
R=CDM_hybrid_output(R,O);
assert(~R.OutputInferenceValid && all(isnan(R.DetailsA(:,3:4)),'all'));
assert(all(isnan(R.DetailsMea1(:,3:4)),'all') && all(isnan(R.DetailsRaw(:,3:4)),'all'));
assert(isfield(R,'DetailsLocalUnchecked') && all(isfinite(R.DetailsLocalUnchecked(:,3))));
assert(all(isnan(R.DetailsLocalUnchecked(:,4))),'Unchecked local SEs acquired p-values.');
% The app's explicit boundary veto must also survive a standalone export.
R.QuadratureStable=true; R.InferenceValid=false;
O.OutputFilePrefix='hybrid_veto'; O.Display=1;
screen=evalc('R=CDM_hybrid_output(R,O);');
assert(~R.OutputInferenceValid,'Explicit inference veto was ignored.');
assert(contains(screen,'Model diagnostics') && ~contains(screen,'Estimation took'),'Saved-fit output printed a spurious estimation time.');
detail='unstable quadrature and stored inference veto suppress public SEs/p-values';
end

function detail=standardCompatibility(outDir)
% Same one-block signature used by ordinary CDM; the second output is optional.
O=struct('Display',0,'OutputDir',outDir,'OutputFilePrefix','ordinary_CDM','xlsOverwrite',1);
R=struct('DetailsA',[1,0,.2,.001;-.4,0,.1,.02],'Dist',[-1;-1],'stats',[]);
N=struct('DetailsA',{{'Constant';'Cost'}}); H=struct('DetailsA',{{'';'tb'}});
head={'Poisson','ordinary CDM compatibility'}; tail=cell(16,2);
[C,path]=genOutput(O,R,head,tail,N,{'DetailsA'},{'DetailsA'},H,{'DetailsA'});
assert(exist(path,'file')==2 && strcmp(path,fullfile(outDir,'ordinary_CDM.xlsx')));
C1=genOutput(O,R,head,tail,N,{'DetailsA'},{'DetailsA'},H,{'DetailsA'});
assert(isequaln(C,C1),'Optional second output changed legacy first output.');
E=readcell(path,'Sheet','Results'); row=find(strcmp(E(:,1),'Cost'),1);
assert(E{row,3}==-.4 && E{row,5}==.1 && strcmp(strtrim(E{row,4}),'**'));
detail='legacy one-output call unchanged; optional filepath is actual saved XLSX';
end

function [R,O]=fixture(outDir)
if ~exist(outDir,'dir'), mkdir(outDir); end
t=(1:24)';
I=struct('Y',1+mod(t,5),'Xa',[ones(24,1),t/24], ...
    'Xstr',(t-mean(t))/std(t),'Xmea',[1+mod(t,3),1+mod(floor(t/2),3),1+mod(floor(t/3),3)]);
E=struct('EvaluateOnly',true,'Display',0,'NB',1,'Truncated',1, ...
    'LatentDistribution','selected','Quadrature','normal-grid','NNodes',64, ...
    'NamesA',{{'Constant';'Travel cost / 100 PLN'}},'NamesMea',{{'Water';'Information';'Sea'}});
R=CDM_hybrid(I,E); p=numel(R.bhat);
if isfield(R,'EstimationSeconds'), R=rmfield(R,'EstimationSeconds'); end
R.bhat(R.Layout.thresholds{1})=[-.8;.4];
R.ihess=diag(linspace(.01,.04,p))+.003*ones(p);
R.std=sqrt(diag(R.ihess));
R.HessDiagnostics.ValidInference=true;
R.QuadratureChecked=true; R.QuadratureStable=true; R.InferenceValid=true;
R.LL=sum(LL_CDM_hybrid(R.bhat,R.INPUT,R.EstimOpt));
logPath=fullfile(outDir,'fixture_log.txt'); fid=fopen(logPath,'w');
fprintf(fid,'Known-covariance test fixture; no estimation performed.\n'); fclose(fid);
O=struct('Display',0,'OutputDir',outDir,'OutputFilePrefix','hybrid_validated', ...
    'xlsOverwrite',1,'SourceScript',[mfilename('fullpath') '.m'],'OutputLogFile',logPath, ...
    'OutputNotes',{{'Travel cost unit','100 PLN';'Xstr scaling','Sample mean zero, sample SD one'}});
end

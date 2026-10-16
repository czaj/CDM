function Results = CDM_demo_hybrid
%CDM_DEMO_HYBRID Reproducible NB2 + ordinal-indicator travel-cost example.
% Requires Optimization Toolbox and Statistics and Machine Learning Toolbox
% (the latter is used only to generate artificial trip counts).
previousRng = rng;
restoreRng = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(20260918,'twister');
n = 1200;
age = randn(n,1);
cost = 0.2 + 2*rand(n,1);
latent = 0.35*age + randn(n,1);
mu = exp(1.4 - 0.8*cost + 0.15*age + 0.45*latent);
alpha = 0.3;
Y = poissrnd(mu .* randg(ones(n,1)/alpha)*alpha);
I = zeros(n,3);
for j = 1:3
    u = rand(n,1);
    propensity = (0.8+0.2*j)*latent + log(u./(1-u));
    I(:,j) = 1 + sum(propensity > [-1 0.5 1.5],2);
end
% Selection is applied AFTER generating both counts and indicators from the
% population model. This is the sampling scheme Truncated=1 corrects for.
keep = Y > 0;
INPUT = struct('Y',Y(keep),'Xa',[ones(sum(keep),1),cost(keep),age(keep)], ...
    'Xstr',age(keep),'Xmea',I(keep,:));
EstimOpt = struct('NB',1,'Truncated',1,'MeaMatrix',ones(3,1), ...
    'NamesA',{{'Constant';'Travel cost';'Age'}},'Display',1);
[INPUT,Results,EstimOpt,OptimOpt] = DataCleanCDM_hybrid(INPUT,EstimOpt);
Results = CDM_hybrid(INPUT,Results,EstimOpt,OptimOpt);
fprintf('Artificial data: true cost coefficient -0.8, fitted %.4f.\n',Results.bhat(2));
fprintf('True latent count coefficient +0.45. Refit with more nodes and check convergence before interpreting.\n');
end

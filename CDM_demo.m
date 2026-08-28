% CDM_demo - a walkthrough of the count data models in this package.
%
% Everything runs on simulated data with known parameters, so you can compare
% the estimates against the values used to generate the data. Run it section by
% section (ctrl+enter) rather than all at once.
%
% Requires the Tools repository (https://github.com/czaj/Tools) on the path.
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

clear
clc
global B_backup

rng(10000001)


%% Simulated recreation demand data
% Xa holds a constant, travel cost and income. The travel cost coefficient is
% negative, so consumer surplus per trip (-1/beta_tc) is positive.

n = 5000; % number of respondents

INPUT.Xa = [ones(n,1),unifrnd(5,50,n,1),unifrnd(-1,1,n,1)]; % Cons, travel cost, income
beta = [1.8;-0.05;0.3];
Lam = exp(INPUT.Xa*beta); % expected number of trips

INPUT.MissingInd = zeros(n,1);

EstimOpt.NP = n;  % number of respondents
EstimOpt.NCT = 1; % observations per respondent
EstimOpt.NamesA = {'Cons';'Travel cost';'Income'};
EstimOpt.TC_par_index = 2; % travel cost is the 2nd parameter of Xa
EstimOpt.ProjectName = 'demo';

fprintf('True beta: %s, implied CS per trip: %.2f\n',mat2str(beta',3),-1/beta(2))


%% Poisson

B_backup = [];
INPUT.Y = poissrnd(Lam);
EstimOpt.NB = 0;
EstimOpt.Zinf = 0;
EstimOpt.Truncated = 0;
EstimOpt.Censored = 0;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.POISS = CDM(INPUT,Results,EstimOpt,OptimOpt);


%% Negative Binomial
% theta = exp(-gamma), so a true theta of 2 means gamma = -log(2) = -0.693.

B_backup = [];
theta = 2;
INPUT.Y = nbinrnd(theta,theta./(theta+Lam));
EstimOpt.NB = 1;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.NB = CDM(INPUT,Results,EstimOpt,OptimOpt);


%% Negative Binomial with covariates of overdispersion
% Overdispersion is allowed to depend on income. A constant is prepended to Xnb
% automatically, so NamesNB only lists the covariates.

B_backup = [];
theta_i = exp(-(0.2 + 0.5*INPUT.Xa(:,3)));
INPUT.Y = nbinrnd(theta_i,theta_i./(theta_i+Lam));
INPUT.Xnb = INPUT.Xa(:,3);
EstimOpt.NB = 1;
EstimOpt.NamesNB = {'Income'};

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.NB_Xnb = CDM(INPUT,Results,EstimOpt,OptimOpt);

INPUT = rmfield(INPUT,'Xnb');
EstimOpt = rmfield(EstimOpt,'NamesNB');


%% Negative Binomial 1
% NB2 lets the variance grow with the square of the mean; NB1 keeps it
% proportional. Which one fits better is an empirical question - CDM_compare
% cannot rank them (they are not nested), so compare AIC or BIC.

B_backup = [];
alpha = 0.8;
theta_nb1 = Lam/alpha; % NB1: theta = lambda/alpha
INPUT.Y = nbinrnd(theta_nb1,theta_nb1./(theta_nb1+Lam));
EstimOpt.NB = 1;
EstimOpt.NBtype = 1;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.NB1 = CDM(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.NBtype = 2; % back to the default for the sections below
EstimOpt.NB = 0;


%% Zero-truncated Poisson (on-site sample)
% Only visitors are observed, so the sample contains no zeros. The draws are
% resampled until positive, which is exactly the truncated distribution.

B_backup = [];
Y = poissrnd(Lam);
while any(Y == 0)
    Y(Y == 0) = poissrnd(Lam(Y == 0));
end
INPUT.Y = Y;
EstimOpt.NB = 0;
EstimOpt.Truncated = 1;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.ZTP = CDM(INPUT,Results,EstimOpt,OptimOpt);


%% Zero-truncated Poisson with endogenous stratification
% On-site sampling where the chance of being intercepted is proportional to the
% number of trips. Under this scheme Y-1 is Poisson distributed.

B_backup = [];
INPUT.Y = poissrnd(Lam) + 1;
EstimOpt.Truncated = 2;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.ZTP_ES = CDM(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.Truncated = 0;


%% Right-censored Poisson
% Counts are only recorded up to a cut-off (e.g. "10 or more trips"). CDM will
% censor the data for you and warn if it had to.

B_backup = [];
INPUT.Y = poissrnd(Lam);
EstimOpt.Censored = 10;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.POISS_cens = CDM(INPUT,Results,EstimOpt,OptimOpt);

EstimOpt.Censored = 0;


%% Zero-inflated Poisson
% A share of the sample never participates, regardless of travel cost. Xzinf
% needs its own constant. (1-q) is the probability of being a potential visitor.

B_backup = [];
% Xzinf holds the covariates only - CDM adds the constant, just as it does for Xnb
INPUT.Xzinf = unifrnd(-1,1,n,1);
INPUT.Xzinf_saved = INPUT.Xzinf; % kept for the model comparison section below
delta = [0.3;-0.8];
q = 1./(1+exp(-[ones(n,1),INPUT.Xzinf]*delta));
INPUT.Y = poissrnd(Lam).*binornd(1,1-q);

EstimOpt.NB = 0;
EstimOpt.Zinf = 1;
EstimOpt.NamesZinf = {'Zvar'};

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.ZIP = CDM(INPUT,Results,EstimOpt,OptimOpt);


%% Zero-inflated Negative Binomial

B_backup = [];
theta = 2;
INPUT.Y = nbinrnd(theta,theta./(theta+Lam)).*binornd(1,1-q);
EstimOpt.NB = 1;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.ZINB = CDM(INPUT,Results,EstimOpt,OptimOpt);

INPUT = rmfield(INPUT,'Xzinf');
EstimOpt.Zinf = 0;


%% Weighted estimation
% Weights are supplied per respondent (or per observation) and are rescaled to
% unit mean. The footer of the output reports weighted maximum likelihood.

B_backup = [];
INPUT.Y = poissrnd(Lam);
INPUT.W = unifrnd(0.5,1.5,n,1);
EstimOpt.NB = 0;

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.POISS_W = CDM(INPUT,Results,EstimOpt,OptimOpt);

INPUT = rmfield(INPUT,'W');


%% Hurdle model
% Two decisions: whether to go at all, and how often. The parts are estimated
% separately because the likelihood separates.

B_backup = [];
delta_h = [0.6;-0.02;0.5];
part = binornd(1,1./(1+exp(-INPUT.Xa*delta_h)));
Ytr = poissrnd(Lam);
while any(Ytr == 0)
    Ytr(Ytr == 0) = poissrnd(Lam(Ytr == 0));
end
INPUT.Y = part.*Ytr;
EstimOpt.NB = 0;
EstimOpt.Truncated = 0;
EstimOpt.NamesH = {'Travel cost';'Income'};

[INPUT,Res0,EstimOpt,OptimOpt] = DataCleanCDM(INPUT,EstimOpt);
Results.CDM0 = Res0.CDM0; % keep the accumulated models, take only the null model
Results.Hurdle = CDM_hurdle(INPUT,Results,EstimOpt,OptimOpt);
Results.Hurdle = CDM_post(Results.Hurdle);


%% Building the travel cost variable
% The money part is straightforward; the value of travel time is a choice, so
% the default returns the whole sensitivity panel.

Trip = struct();
Trip.Distance = unifrnd(5,200,n,1);   % one way, km
Trip.Time = Trip.Distance/70;         % one way, hours
Trip.CostPerKm = 0.35;
Trip.Wage = unifrnd(15,60,n,1);
Trip.ValuesTime = binornd(1,0.6,n,1); % 1 = would have preferred a shorter trip

[TCpanel,TCdetail] = TravelCost(Trip);                              % fractions of the wage
TCfull = TravelCost(Trip,struct('TimeRule','full'));                % full wage
TCmix  = TravelCost(Trip,struct('TimeRule','full_or_zero'));        % full wage, or nothing


%% Post-estimation: consumer surplus, marginal effects, fit
% Everything a travel cost study reports, with Krinsky-Robb intervals.

Results.POISS = CDM_post(Results.POISS);


%% Comparing models against each other
% Specification tests only make sense on one and the same sample, so this section
% fits the competing models to a single dataset - here one with genuine
% overdispersion and an excess of zeros.

B_backup = [];
thc = 1.5;
qc = 1./(1+exp(-[ones(n,1),INPUT.Xzinf_saved]*[0.3;-0.8]));
Yc = nbinrnd(thc,thc./(thc+Lam)).*binornd(1,1-qc);

INPUTc = struct('Xa',INPUT.Xa,'Y',Yc,'MissingInd',zeros(n,1));
EstimOptC = EstimOpt;
EstimOptC.NB = 0; EstimOptC.Zinf = 0;
[INc,Res0,EstimOptC,OptimOptC] = DataCleanCDM(INPUTc,EstimOptC);
Cmp.CDM0 = Res0.CDM0;
Cmp.POISS = CDM(INc,Cmp,EstimOptC,OptimOptC);

B_backup = [];
EstimOptC.NB = 1;
[INc,Res0,EstimOptC,OptimOptC] = DataCleanCDM(INPUTc,EstimOptC);
Cmp.CDM0 = Res0.CDM0;
Cmp.NB = CDM(INc,Cmp,EstimOptC,OptimOptC);

B_backup = [];
INPUTc.Xzinf = INPUT.Xzinf_saved;
EstimOptC.NB = 0; EstimOptC.Zinf = 1; EstimOptC.NamesZinf = {'Zvar'};
[INc,Res0,EstimOptC,OptimOptC] = DataCleanCDM(INPUTc,EstimOptC);
Cmp.CDM0 = Res0.CDM0;
Cmp.ZIP = CDM(INc,Cmp,EstimOptC,OptimOptC);

CDM_compare(Cmp.POISS,Cmp.NB);   % boundary-corrected LR: is overdispersion real?
CDM_compare(Cmp.POISS,Cmp.ZIP);  % Voung: is zero inflation needed?


%% Comparing models
% stats = [LL; LL0; McFadden R2; (unused); AIC/n; BIC/n; n; r; k]

fprintf('\n%-28s %12s %10s %10s\n','model','LL','AIC/n','k')
models = {'POISS','NB','NB1','ZIP','ZINB','Hurdle'};
for i = 1:numel(models)
    if isfield(Results,models{i})
        s = Results.(models{i}).stats;
        fprintf('%-28s %12.2f %10.4f %10d\n',models{i},s(1),s(5),s(9))
    end
end

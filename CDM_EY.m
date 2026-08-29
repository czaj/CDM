function EY = CDM_EY(b,Xa,Xnb,Xzinf,EstimOpt)
% CDM_EY expected number of trips implied by a fitted count data model.
%
% Syntax:   EY = CDM_EY(b,Xa,Xnb,Xzinf,EstimOpt)
%
% The expectation matches the sampling scheme the model was estimated under, which
% is what makes it usable for welfare measures:
%
%   Truncated = 0   E[Y]      = lambda                  general population sample
%   Truncated = 1   E[Y|Y>0]  = lambda/(1-P(0))         on-site sample
%   Truncated = 2   E*[Y]     = E[Y^2]/E[Y]             on-site, endogenous stratification
%
% Zero inflation scales the expectation by the probability of being a participant.
% Right censoring describes how trips were recorded rather than how many were
% taken, so it does not enter here - welfare measures need the latent expectation.
%
% See also CDM, CDM_post
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

NVarA = EstimOpt.NVarA;
Lam = exp(Xa*b(1:NVarA));

if EstimOpt.NB == 1
    % theta = lambda^(2-P)/alpha: P = 2 is NB2, P = 1 is NB1, P free is NB-P
    if isfield(EstimOpt,'NBtype') && EstimOpt.NBtype == 3
        powerTheta = 2 - b(NVarA+EstimOpt.NVarNB+1);
    elseif isfield(EstimOpt,'NBtype') && EstimOpt.NBtype == 1
        powerTheta = 1;
    else
        powerTheta = 0;
    end
    theta = exp(powerTheta*log(Lam) - Xnb*b(NVarA+1:NVarA+EstimOpt.NVarNB));
else
    theta = [];
end

switch EstimOpt.Truncated
    case 0
        EY = Lam;
    case 1 % zero-truncated: divide by the probability of at least one trip
        if isempty(theta)
            oneMinusP0 = -expm1(-Lam);
        else
            oneMinusP0 = -expm1(theta.*log(theta./(theta+Lam)));
        end
        EY = Lam./oneMinusP0;
    case 2 % size-biased: E*[Y] = E[Y^2]/E[Y] = 1 + lambda (+ overdispersion)
        if isempty(theta)
            EY = 1 + Lam;
        else
            EY = 1 + Lam + Lam./theta;
        end
    otherwise
        error('Unsupported EstimOpt.Truncated')
end

if EstimOpt.Zinf == 1
    bzinf = b(NVarA+EstimOpt.NVarNB+1:NVarA+EstimOpt.NVarNB+EstimOpt.NVarZinf);
    if isfield(EstimOpt,'probit') && EstimOpt.probit == 1
        q = normcdf(Xzinf*bzinf);
    else
        q = 1./(1+exp(-Xzinf*bzinf));
    end
    EY = (1-q).*EY;
end
end

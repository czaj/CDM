function [f,g] = LL_CDM(Y,Xa,Xnb, Xzinf,EstimOpt,b)
%upper censoring possible -> Censored = upper limit
%possible truncation Truncated = 1
%possible endogenous truncation (trunc + strat) Truncated = 2
%WIP! Zinf = 1 makes the model zero inflated (logit is default, probit in
%progress (probit ==1)

Zinf = EstimOpt.Zinf;
NB = EstimOpt.NB;
Truncated = EstimOpt.Truncated;
NVarA = EstimOpt.NVarA;
NVarNB = EstimOpt.NVarNB;
NVarZinf = EstimOpt.NVarZinf;
Censored = EstimOpt.Censored;
if isfield(EstimOpt,'NBtype') && ~isempty(EstimOpt.NBtype)
    NBtype = EstimOpt.NBtype; % 2 = NB2 (default), 1 = NB1, 3 = NB-P
else
    NBtype = 2;
end

% All three parameterisations are NB2 with theta = lambda^(2-P)/alpha:
%   NB2   P = 2   Var = lambda + alpha*lambda^2
%   NB1   P = 1   Var = lambda*(1+alpha)
%   NB-P  P free  Var = lambda + alpha*lambda^P
% powerTheta is 2-P, the exponent with which lambda enters theta. It is the only
% thing that distinguishes them, in the likelihood and in the gradient alike.
nExtra = 0;
powerTheta = 0;
if NB == 1
    switch NBtype
        case 1
            powerTheta = 1;
        case 3
            nExtra = 1; % P is estimated, and sits right after the alpha block
            powerTheta = 2 - b(NVarA+NVarNB+1);
        otherwise
            powerTheta = 0;
    end
end
RealMin = EstimOpt.RealMin; 

% Censoring masks, built once - they were re-evaluated on every use before.
if Censored > 0
    iC = (Y == Censored); % observations at the censoring point
    iU = ~iC;             % observations below it
else
    iC = false(size(Y));
    iU = true(size(Y));
end

%loglik
if NB == 0
    if Censored ==0
        fit = Xa*b(1:NVarA);
        if Truncated == 0
            f = (fit).*Y - gammaln(Y+1) - exp(fit);
        elseif Truncated == 2
            f = (fit).*(Y-1) - gammaln(Y) - exp(fit);
        elseif Truncated ==1
            f = (fit).*Y - gammaln(Y+1) - exp(fit) - log(-expm1(-exp(fit)));
        end
    else % Censored >0
        f = zeros(length(Y),1);

        ba = b(1:NVarA);
        cXa = Xa(iC,:);
        fitc = cXa*ba;   % linear index at the censoring point
        clam = exp(fitc);
        fitu = Xa(iU,:)*ba; % linear index below it
        lamu = exp(fitu);
        Yu = Y(iU);
        
        if Truncated == 0
            %if Y< C
            f(iU) = fitu.*Yu - gammaln(Yu+1) - lamu;

            %if Y=C
            vec = 0:(Censored - 1);
            vec = vec(ones(length(clam),1),:);
            fitc = fitc(:,ones(Censored,1));
            ftmp = gammainc(clam,Censored,'lower'); % = 1 - sum_{k<C} Poisson(k)

            if RealMin == 0
                f(iC) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp);
            end

        elseif Truncated == 2
            %if Y< C
            f(iU) = fitu.*(Yu-1) - gammaln(Yu) - lamu;

            %if Y=C
            vec = 0:(Censored - 2);
            vec = vec(ones(length(clam),1),:);

            fitc = fitc(:,ones(Censored-1,1));
            ftmp = gammainc(clam,Censored-1,'lower'); % Y-1 is Poisson under stratification

            if RealMin == 0
                f(iC) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp);
            end

        elseif Truncated ==1
            f(iU) = fitu.*Yu - gammaln(Yu+1) - lamu - log(-expm1(-lamu));
            
            %for Y=C
            vec = 0:(Censored - 1);
            vec = vec(ones(length(clam),1),:);
            fitc = fitc(:,ones(Censored,1));
            ftmp = gammainc(clam,Censored,'lower');

            if RealMin == 0
                f(iC) = log(ftmp) - log(-expm1(-clam));
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp) - log(-expm1(-clam));
            end
            
        end

    end


elseif NB == 1
    fit = Xa*b(1:NVarA);
    % fit2 is -log(theta) throughout, which is what lets all three
    % parameterisations share every branch below.
    fit2 = Xnb*b(NVarA+1:NVarA+NVarNB) - powerTheta*fit;
    theta = exp(-fit2);

    if Censored ==0
        if Truncated == 0
            f = gammaln(theta+Y)-gammaln(theta)-gammaln(Y+1)-theta.*fit2 + Y.*fit - (Y+theta).*log(theta + exp(fit));

        elseif Truncated == 1
            ThetaLam = theta./(exp(fit) + theta);

            f = gammaln(theta+Y)-gammaln(theta)-gammaln(Y+1)-theta.*fit2 + Y.*fit - (Y+theta).*log(theta + exp(fit)) - log(-expm1(theta.*log(ThetaLam)));
        elseif Truncated == 2

            f = gammaln(theta+Y)-gammaln(theta)-gammaln(Y+1)-theta.*fit2 + Y.*fit - (Y+theta).*log(theta + exp(fit)) + log(Y) - fit;
        end

    else % Censored >0
        f = zeros(length(Y),1);

        cfit = fit(iC);
        cfit2 = fit2(iC);
        ctheta = theta(iC);

        vec = 0:(Censored - 1);
        vec = vec(ones(length(cfit),1),:);

        % P(Y >= C) for a Negative Binomial is the regularized incomplete beta.
        % Computing it directly avoids the cancellation in 1 - sum_k P(k), which
        % underflows to exactly zero once the censoring point sits far in the tail
        % and then feeds -Inf to the optimizer.
        cLam1 = exp(cfit);
        cTh1 = ctheta;
        survC = betainc(cLam1./(cTh1+cLam1),Censored,cTh1);

        cfit = cfit(:,ones(Censored,1));
        cfit2 = cfit2(:,ones(Censored,1));
        ctheta = ctheta(:,ones(Censored,1));


        if Truncated == 0
            %if Y< C
            f(iU) = gammaln(theta(iU)+Y(iU))-gammaln(theta(iU))-gammaln(Y(iU)+1)-theta(iU).*fit2(iU) ...
                + Y(iU).*fit(iU) - (Y(iU)+theta(iU)).*log(theta(iU) + exp(fit(iU)));

            %if Y= C

            ftmp = survC;

            if RealMin == 0
                f(iC) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp);
            end

        elseif Truncated == 2
            %if Y< C
            f(iU) = gammaln(theta(iU)+Y(iU))-gammaln(theta(iU))-gammaln(Y(iU)+1)...
                -theta(iU).*fit2(iU) + Y(iU).*fit(iU)...
                - (Y(iU)+theta(iU)).*log(theta(iU) + exp(fit(iU)))...
                + log(Y(iU)) - fit(iU);

            %if Y= C
            ftmp = 1 - sum(exp(gammaln(ctheta+vec)- gammaln(ctheta) - gammaln(vec+1)...
                - ctheta.*cfit2 + vec .* cfit - (ctheta+vec).*log(ctheta+exp(cfit)) + log(vec) - cfit),2);

            if RealMin == 0
                f(iC) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp);
            end


        elseif Truncated == 1
            %if Y< C
%             cfit = cfit(:,1:end-1);
%             cfit2 = cfit2(:,1:end-1);
%             ctheta = ctheta(:,1:end-1);

            ThetaLam = theta./(exp(fit) + theta);
            f(iU) = gammaln(theta(iU)+Y(iU))-gammaln(theta(iU))-gammaln(Y(iU)+1)...
                -theta(iU).*fit2(iU) + Y(iU).*fit(iU)...
                - (Y(iU)+theta(iU)).*log(theta(iU) + exp(fit(iU)))...
                - log(-expm1(theta(iU).*log(ThetaLam(iU))));
            
            %if Y=C
            % - log(1-(ThetaLam).^(theta))
            
%             cThetaLam = ThetaLam(iC);
            trunc_corr = -expm1(theta(iC).*log(ThetaLam(iC)));
            ftmp = survC;
           
%                 vec = 1:(Censored - 1);
%             vec = vec(ones(length(cfit),1),:);
%             ftmp = 1 - sum(exp(gammaln(ctheta+vec)-gammaln(ctheta)-gammaln(vec+1)...
%                 -ctheta.*cfit2 + vec.*cfit - (vec+ctheta).*log(ctheta + exp(cfit))),2)./trunc_corr;
    

% Without censoring
%             ThetaLam = theta./(exp(fit) + theta);
%             f = gammaln(theta+Y)-gammaln(theta)-gammaln(Y+1)-theta.*fit2 + Y.*fit - (Y+theta).*log(theta + exp(fit)) - log(-expm1(theta.*log(ThetaLam)));

% Without truncation
%             f(iU) = gammaln(theta(iU)+Y(iU))-gammaln(theta(iU))-gammaln(Y(iU)+1)-theta(iU).*fit2(iU) ...
%                 + Y(iU).*fit(iU) - (Y(iU)+theta(iU)).*log(theta(iU) + exp(fit(iU)));
% 
%             ftmp = 1 - sum( exp( gammaln(ctheta+vec)- gammaln(ctheta) - gammaln(vec+1)...
%                 - ctheta.*cfit2 + vec .* cfit - (ctheta+vec).*log(ctheta+exp(cfit))),2);

            if RealMin == 0
                f(iC) = log(ftmp) - log(trunc_corr);
            else
                ftmp = max(ftmp,realmin);
                f(iC) = log(ftmp)  - log(trunc_corr);
            end

        end
% elseif NB==2 not coded, but is in Wiktor's code
    end
    
end

if Zinf == 1
    bzinf = b(NVarA+NVarNB+nExtra+1:NVarA+NVarNB+nExtra+NVarZinf);
    if isfield(EstimOpt, 'probit') && EstimOpt.probit == 1
        q = normcdf(Xzinf*bzinf);
    else
        q = exp(Xzinf*bzinf)./(1+exp(Xzinf*bzinf));
    end
    f0 = f(Y == 0);
    f(Y == 0) = log((1-q(Y==0)).*exp(f0) + q(Y==0));
    f(Y > 0) = log(1-q(Y>0)) + f(Y>0);

end
f = -f;


if nargout == 2
    Lam = exp(Xa*b(1:NVarA));
    if NB == 0
        if Censored ==0
            if Truncated ==0
                g = Xa.*Y - Lam.*Xa;
                g = -g;
            elseif  Truncated ==1
                truncation_correction = Lam .* exp(-Lam) ./ (-expm1(-Lam));
                g = Xa .* (Y - Lam - truncation_correction);
                g = -g;
%                 f = (fit).*Y - gammaln(Y+1) - exp(fit) - log(-expm1(-exp(fit)));
            elseif  Truncated ==2
                g = Xa.*(Y - 1) - Lam.*Xa;
                g = -g;
            end
        else % Censored >0
            g = zeros(length(Y),EstimOpt.NVarA);

            if Truncated ==0

                %for Y<C
                g(iU,:) = Xa(iU,:).*Y(iU) - Lam(iU).*Xa(iU,:);
                %for Y=C
                cLam = Lam(iC,:);
                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored,1))-vec),2);
                gtmp = gtmp./ftmp;

                g(iC,:) = gtmp.*cXa;

                g = -g;

            elseif  Truncated ==2
                %for Y<C

                g(iU,:) = Xa(iU,:).*(Y(iU)-1) - Lam(iU).*Xa(iU,:);
                %for Y=C
                cLam = Lam(iC,:);

                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored-1,1))-vec),2);
                gtmp = gtmp./ftmp;

                g(iC,:) = gtmp.*cXa;
                g = -g;


            elseif  Truncated ==1
                truncation_correction = Lam .* exp(-Lam) ./ (-expm1(-Lam));

                %for Y<C
                g(iU,:) = Xa(iU,:) .* (Y(iU,:) - Lam(iU,:) - truncation_correction(iU,:));
                
                %for Y=C
                cLam = Lam(iC,:);
            
                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .*(cLam(:,ones(Censored,1))-vec),2);
%                 gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored,1))-vec),2);
%                 ftmp = (1 - sum(exp(vec.*fitc - clam - gammaln(vec+1)),2))./(1 - exp(-clam));
                gtmp = gtmp./ftmp  - truncation_correction(iC,:);

                g(iC,:) = gtmp.*cXa;

                g = -g;

            end

        end
        
    else %if NB==1
        if Censored ==0
%         if Truncated == 0
            %derivative with respect to beta
            g = Xa.*Y - (Y + theta) .* Lam.*Xa./(theta+Lam) ;
            % dnb collects everything that reaches the likelihood through theta.
            % The gamma-gradient is dnb.*Xnb because d(log theta)/d(gamma) = -Xnb;
            % under NB1 the same channel is open to beta, where d(log theta)/d(beta) = Xa.
            dnb  = theta.*(-psi(theta+Y) + psi(theta)-1 + fit2 + log(theta+Lam) + (theta+Y)./(theta+Lam));

            if Truncated == 2
                g = g - Xa;

            elseif Truncated ==1

                gtmp = (ThetaLam.^(theta+1))./(-expm1(theta.*log(ThetaLam)));
                gtmp = gtmp.*Xa.*Lam;
                g = g - gtmp;

                dtr = (ThetaLam.^theta).*(Lam + (theta+Lam).*log(ThetaLam));
                dtr = dtr./((theta+Lam).*(-expm1(theta.*log(ThetaLam))));
                dnb = dnb - dtr.*theta;
            end
            dTheta = dnb; % gnb = dTheta.*Xnb, since d(log theta)/d(gamma) = -Xnb

        else % Censored >0
            txnb = theta.*Xnb;

            % NB pmf at the censored cell, for k = 0 ... Censored-1. Everything goes
            % through gammaln: gamma(theta) overflows to Inf for theta > 171, which is
            % reached at ln(alpha) < -5.2, i.e. exactly when the data are near-Poisson.
            clam = exp(cfit);
            u = ctheta./(ctheta+clam);
            fk = exp(gammaln(ctheta+vec) - gammaln(ctheta) - gammaln(vec+1) ...
                + ctheta.*log(u) + vec.*cfit - vec.*log(ctheta+clam));

            % d log f(k) / d beta = u*(k-lambda)                  (times Xa)
            % d log f(k) / d s    = -theta*dlnf_ds, s = Xnb*gamma (times Xnb)
            dlnf_ds = psi(ctheta+vec) - psi(ctheta) + log(u) + 1 - u - vec./(ctheta+clam);

            % Y<C
            g = zeros(length(Y),EstimOpt.NVarA);
            g(iU,:) = Xa(iU,:).*Y(iU) - (Y(iU)...
                + theta(iU)).*Lam(iU).*Xa(iU,:)./(theta(iU)+Lam(iU)) ;
            
            gnb = zeros(length(Y),EstimOpt.NVarNB);
            gnb_tmp  = -psi(theta(iU)+Y(iU))+psi(theta(iU))-1+fit2(iU)...
                +log(theta(iU)+Lam(iU))+ (theta(iU)+Y(iU))./(theta(iU)+Lam(iU));
            dTheta = zeros(length(Y),1);
            dTheta(iU) = theta(iU).*gnb_tmp;
                
            if Truncated == 2
                % Endogenous stratification (on-site sample): the censored cell mixes
                % the size-biased masses P*(k) = k*f(k)/lambda, not f(k). Both the
                % weights and the per-k derivative w.r.t. beta change accordingly:
                %   d log P*(k) / d beta = u*(k-lambda) - 1
                %   d log P*(k) / d s    = d log f(k) / d s   (log k - log lambda is free of gamma)
                pk = fk.*vec./clam;

                % Y<C
                g(iU,:) = g(iU,:) - Xa(iU,:);

                % Y=C: d log ftmp / d beta = [ sum_k P*(k) + u*sum_k P*(k)*(lambda-k) ] / ftmp
                gtmp = (sum(pk,2) + u(:,1).*sum(pk.*(clam-vec),2))./ftmp;
                g(iC,:) = gtmp.*Xa(iC,:);

                % Y=C: d log ftmp / d s = theta * sum_k P*(k)*dlnf_ds / ftmp
                gtmp = ctheta(:,1).*sum(pk.*dlnf_ds,2)./ftmp;
                dTheta(iC) = gtmp;

            else % Truncated == 0 or 1 - the censored cell mixes the plain NB masses
                % d log ftmp / d beta = u * sum_k f(k)*(lambda-k) / ftmp
                gtmp = u(:,1).*sum(fk.*(clam-vec),2)./ftmp;
                g(iC,:) = gtmp.*Xa(iC,:);

                % d log ftmp / d s = theta * sum_k f(k)*dlnf_ds / ftmp
                gtmp = ctheta(:,1).*sum(fk.*dlnf_ds,2)./ftmp;
                dTheta(iC) = gtmp;

                if Truncated == 1
                    % The zero-truncation correction -log(1-u^theta) has the same form
                    % for both groups, so it is subtracted from every row.
                    gtrunc = (ThetaLam.^(theta+1))./(-expm1(theta.*log(ThetaLam)));
                    gtrunc = gtrunc.*Xa.*Lam;
                    g = g - gtrunc;

                    gtruncnb = (ThetaLam.^theta).*(Lam + (theta+Lam).*log(ThetaLam));
                    gtruncnb = gtruncnb./((theta+Lam).*(-expm1(theta.*log(ThetaLam))));
                    dTheta = dTheta - gtruncnb.*theta;
                end
            end
        end


    gnb = dTheta.*Xnb;
    if powerTheta ~= 0 % lambda enters theta, so beta reaches it too
        g = g - powerTheta*dTheta.*Xa;
    end
    if nExtra == 1 % d(log theta)/d(P) = -log(lambda)
        gP = dTheta.*fit;
    else
        gP = zeros(length(Y),0);
    end

    g = -[g, gnb, gP];
%     g

    end
    
    if Zinf ==1
    %for probit
        if isfield(EstimOpt, 'probit') && EstimOpt.probit == 1
            gzinf = normpdf(Xzinf*bzinf).*Xzinf;
        else
            gzinf = q.*(1-q).*Xzinf;
        end
        L = exp(-f);
        g(Y == 0,:) = exp(f0).*(1-q(Y == 0)).*g(Y == 0,:)./L(Y == 0);
        gzinf(Y>0,:) = gzinf(Y>0,:)./(1-q(Y > 0));
        gzinf(Y == 0,:) = gzinf(Y == 0,:).*(exp(f0)-1)./L(Y == 0);
        g = [g, gzinf];
    end
    
end



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
RealMin = EstimOpt.RealMin; 

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

        cXa = Xa(Y == Censored,:);
        clam = exp(cXa*b(1:NVarA));

        fitc = cXa*b(1:NVarA);
        
        if Truncated == 0
            %if Y< C
            f(Y ~= Censored) = (Xa(Y ~= Censored,:)*b(1:NVarA)).*Y(Y ~= Censored) - gammaln(Y(Y ~= Censored)+1) - exp(Xa(Y ~= Censored,:)*b(1:NVarA));

            %if Y=C
            vec = 0:(Censored - 1);
            vec = vec(ones(length(clam),1),:);
            fitc = fitc(:,ones(Censored,1));
            ftmp = 1 - sum(exp(vec.*fitc - clam - gammaln(vec+1)),2);

            if RealMin == 0
                f(Y == Censored) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp);
            end

        elseif Truncated == 2
            %if Y< C
            f(Y ~= Censored) = (Xa(Y ~= Censored,:)*b(1:NVarA)).*(Y(Y ~= Censored)-1) - gammaln(Y(Y ~= Censored)) - exp(Xa(Y ~= Censored,:)*b(1:NVarA));

            %if Y=C
            vec = 0:(Censored - 2);
            vec = vec(ones(length(clam),1),:);

            fitc = cXa*b(1:NVarA);
            fitc = fitc(:,ones(Censored-1,1));
            ftmp = 1 - sum(exp( vec.*fitc - clam - gammaln(vec+1)),2);

            if RealMin == 0
                f(Y == Censored) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp);
            end

        elseif Truncated ==1
            f(Y ~= Censored) = (Xa(Y ~= Censored,:)*b(1:NVarA)).*Y(Y ~= Censored) - gammaln(Y(Y ~= Censored)+1) - exp(Xa(Y ~= Censored,:)*b(1:NVarA)) - log(-expm1(-exp(Xa(Y ~= Censored,:)*b(1:NVarA))));
            
            %for Y=C
            vec = 0:(Censored - 1);
            vec = vec(ones(length(clam),1),:);
            fitc = fitc(:,ones(Censored,1));
            ftmp = (1 - sum(exp(vec.*fitc - clam - gammaln(vec+1)),2));

            if RealMin == 0
                f(Y == Censored) = log(ftmp) - log(-expm1(-clam));
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp) - log(-expm1(-clam));
            end
            
        end

    end


elseif NB == 1
    theta = exp(-Xnb*b(NVarA+1:NVarA+NVarNB));
    fit = Xa*b(1:NVarA);
    fit2 = Xnb*b(NVarA+1:NVarA+NVarNB);

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

        cfit = fit(Y == Censored);
        cfit2 = fit2(Y == Censored);
        ctheta = theta(Y == Censored);

        vec = 0:(Censored - 1);
        vec = vec(ones(length(cfit),1),:);

        cfit = cfit(:,ones(Censored,1));
        cfit2 = cfit2(:,ones(Censored,1));
        ctheta = ctheta(:,ones(Censored,1));


        if Truncated == 0
            %if Y< C
            f(Y ~= Censored) = gammaln(theta(Y ~= Censored)+Y(Y ~= Censored))-gammaln(theta(Y ~= Censored))-gammaln(Y(Y ~= Censored)+1)-theta(Y ~= Censored).*fit2(Y ~= Censored) ...
                + Y(Y ~= Censored).*fit(Y ~= Censored) - (Y(Y ~= Censored)+theta(Y ~= Censored)).*log(theta(Y ~= Censored) + exp(fit(Y ~= Censored)));

            %if Y= C

            ftmp = 1 - sum( exp( gammaln(ctheta+vec)- gammaln(ctheta) - gammaln(vec+1)...
                - ctheta.*cfit2 + vec .* cfit - (ctheta+vec).*log(ctheta+exp(cfit))),2);

            if RealMin == 0
                f(Y == Censored) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp);
            end

        elseif Truncated == 2
            %if Y< C
            f(Y ~= Censored) = gammaln(theta(Y ~= Censored)+Y(Y ~= Censored))-gammaln(theta(Y ~= Censored))-gammaln(Y(Y ~= Censored)+1)...
                -theta(Y ~= Censored).*fit2(Y ~= Censored) + Y(Y ~= Censored).*fit(Y ~= Censored)...
                - (Y(Y ~= Censored)+theta(Y ~= Censored)).*log(theta(Y ~= Censored) + exp(fit(Y ~= Censored)))...
                + log(Y(Y ~= Censored)) - fit(Y ~= Censored);

            %if Y= C
            ftmp = 1 - sum(exp(gammaln(ctheta+vec)- gammaln(ctheta) - gammaln(vec+1)...
                - ctheta.*cfit2 + vec .* cfit - (ctheta+vec).*log(ctheta+exp(cfit)) + log(vec) - cfit),2);

            if RealMin == 0
                f(Y == Censored) = log(ftmp);
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp);
            end


        elseif Truncated == 1
            %if Y< C
%             cfit = cfit(:,1:end-1);
%             cfit2 = cfit2(:,1:end-1);
%             ctheta = ctheta(:,1:end-1);

            ThetaLam = theta./(exp(fit) + theta);
            f(Y ~= Censored) = gammaln(theta(Y ~= Censored)+Y(Y ~= Censored))-gammaln(theta(Y ~= Censored))-gammaln(Y(Y ~= Censored)+1)...
                -theta(Y ~= Censored).*fit2(Y ~= Censored) + Y(Y ~= Censored).*fit(Y ~= Censored)...
                - (Y(Y ~= Censored)+theta(Y ~= Censored)).*log(theta(Y ~= Censored) + exp(fit(Y ~= Censored)))...
                - log(-expm1(theta(Y ~= Censored).*log(ThetaLam(Y ~= Censored))));
            
            %if Y=C
            % - log(1-(ThetaLam).^(theta))
            
%             cThetaLam = ThetaLam(Y == Censored);
            trunc_corr = -expm1(theta(Y == Censored).*log(ThetaLam(Y == Censored)));
            ftmp = 1 - sum(exp(gammaln(ctheta+vec)-gammaln(ctheta)-gammaln(vec+1)...
                -ctheta.*cfit2 + vec.*cfit - (vec+ctheta).*log(ctheta + exp(cfit))),2);
           
%                 vec = 1:(Censored - 1);
%             vec = vec(ones(length(cfit),1),:);
%             ftmp = 1 - sum(exp(gammaln(ctheta+vec)-gammaln(ctheta)-gammaln(vec+1)...
%                 -ctheta.*cfit2 + vec.*cfit - (vec+ctheta).*log(ctheta + exp(cfit))),2)./trunc_corr;
    

% Without censoring
%             ThetaLam = theta./(exp(fit) + theta);
%             f = gammaln(theta+Y)-gammaln(theta)-gammaln(Y+1)-theta.*fit2 + Y.*fit - (Y+theta).*log(theta + exp(fit)) - log(-expm1(theta.*log(ThetaLam)));

% Without truncation
%             f(Y ~= Censored) = gammaln(theta(Y ~= Censored)+Y(Y ~= Censored))-gammaln(theta(Y ~= Censored))-gammaln(Y(Y ~= Censored)+1)-theta(Y ~= Censored).*fit2(Y ~= Censored) ...
%                 + Y(Y ~= Censored).*fit(Y ~= Censored) - (Y(Y ~= Censored)+theta(Y ~= Censored)).*log(theta(Y ~= Censored) + exp(fit(Y ~= Censored)));
% 
%             ftmp = 1 - sum( exp( gammaln(ctheta+vec)- gammaln(ctheta) - gammaln(vec+1)...
%                 - ctheta.*cfit2 + vec .* cfit - (ctheta+vec).*log(ctheta+exp(cfit))),2);

            if RealMin == 0
                f(Y == Censored) = log(ftmp) - log(trunc_corr);
            else
                ftmp = max(ftmp,realmin);
                f(Y == Censored) = log(ftmp)  - log(trunc_corr);
            end

        end
% elseif NB==2 not coded, but is in Wiktor's code
    end
    
end

if Zinf == 1
    bzinf = b(NVarA+1+NVarNB:NVarA + NVarNB + NVarZinf);
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
                g = Xa.*Y(:, ones(1, NVarA)) - Lam(:, ones(1, NVarA)).*Xa;
                g = -g;
            elseif  Truncated ==1
                truncation_correction = Lam .* exp(-Lam) ./ (-expm1(-Lam));
                g = Xa .* (Y - Lam - truncation_correction);
                g = -g;
%                 f = (fit).*Y - gammaln(Y+1) - exp(fit) - log(-expm1(-exp(fit)));
            elseif  Truncated ==2
                g = Xa.*(Y(:, ones(1, NVarA)) - 1) - Lam(:, ones(1, NVarA)).*Xa;
                g = -g;
            end
        else % Censored >0
            g = zeros(length(Y),EstimOpt.NVarA);

            if Truncated ==0

                %for Y<C
                g(Y ~= Censored,:) = Xa(Y ~= Censored,:).*Y(Y ~= Censored, ones(1, NVarA)) - Lam(Y ~= Censored, ones(1, NVarA)).*Xa(Y ~= Censored,:);
                %for Y=C
                cLam = Lam(Y == Censored,:);
                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored,1))-vec),2);
                gtmp = gtmp./ftmp;

                g(Y == Censored,:) = gtmp(:, ones(1, NVarA)).*cXa;

                g = -g;

            elseif  Truncated ==2
                %for Y<C

                g(Y ~= Censored,:) = Xa(Y ~= Censored,:).*(Y(Y ~= Censored, ones(1, NVarA))-1) - Lam(Y ~= Censored, ones(1, NVarA)).*Xa(Y ~= Censored,:);
                %for Y=C
                cLam = Lam(Y == Censored,:);

                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored-1,1))-vec),2);
                gtmp = gtmp./ftmp;

                g(Y == Censored,:) = gtmp(:, ones(1, NVarA)).*cXa;
                g = -g;


            elseif  Truncated ==1
                truncation_correction = Lam .* exp(-Lam) ./ (-expm1(-Lam));

                %for Y<C
                g(Y ~= Censored,:) = Xa(Y ~= Censored,:) .* (Y(Y ~= Censored,:) - Lam(Y ~= Censored,:) - truncation_correction(Y ~= Censored,:));
                
                %for Y=C
                cLam = Lam(Y == Censored,:);
            
                gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .*(cLam(:,ones(Censored,1))-vec),2);
%                 gtmp = sum(exp(vec.*fitc -cLam-gammaln(vec+1)) .* (cLam(:,ones(Censored,1))-vec),2);
%                 ftmp = (1 - sum(exp(vec.*fitc - clam - gammaln(vec+1)),2))./(1 - exp(-clam));
                gtmp = gtmp./ftmp  - truncation_correction(Y == Censored,:);

                g(Y == Censored,:) = gtmp(:, ones(1, NVarA)).*cXa;

                g = -g;

            end

        end
        
    else %if NB==1
        if Censored ==0
%         if Truncated == 0
            %derivative with respect to beta
            g = Xa.*Y(:, ones(1, NVarA)) - (Y(:, ones(1, NVarA)) + theta(:, ones(1, NVarA))) .* Lam(:, ones(1, NVarA)).*Xa./(theta(:, ones(1, NVarA))+Lam(:, ones(1, NVarA))) ;
            %derivative with respect to beta nb
            txnb = theta(:, ones(1, NVarNB)).*Xnb; % d(theta)/d(b)*(-1)
            fit2 = Xnb*b(NVarA+1:NVarA+NVarNB); 
            gnb  = -psi(theta+Y) + psi(theta)-1 + fit2 + log(theta+Lam) + (theta+Y)./(theta+Lam); %gnb*(-1)/txnb
            gnb  = txnb.*gnb(:, ones(1, NVarNB));

            if Truncated == 2
                g = g - Xa;
    
            elseif Truncated ==1
       
                gtmp = (ThetaLam.^(theta+1))./(-expm1(theta.*log(ThetaLam)));
                gtmp = gtmp.*Xa.*Lam;
                g = g - gtmp;
    
                gtmp = (ThetaLam.^theta).*(Lam + (theta+Lam).*log(ThetaLam));
                gtmp = gtmp./((theta+Lam).*(-expm1(theta.*log(ThetaLam))));
                gtmp = gtmp.*Xnb.*theta;
                gnb = gnb - gtmp;
            end

        else % Censored >0
            txnb = theta(:, ones(1, NVarNB)).*Xnb;
            fit2 = Xnb*b(NVarA+1:NVarA+NVarNB);

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
            g(Y ~= Censored,:) = Xa(Y ~= Censored,:).*Y(Y ~= Censored, ones(1, NVarA)) - (Y(Y ~= Censored, ones(1, NVarA))...
                + theta(Y ~= Censored, ones(1, NVarA))).*Lam(Y ~= Censored, ones(1, NVarA)).*Xa(Y ~= Censored,:)./(theta(Y ~= Censored, ones(1, NVarA))+Lam(Y ~= Censored, ones(1, NVarA))) ;
            
            gnb = zeros(length(Y),EstimOpt.NVarNB);
            gnb_tmp  = -psi(theta(Y ~= Censored)+Y(Y ~= Censored))+psi(theta(Y ~= Censored))-1+fit2(Y ~= Censored)...
                +log(theta(Y ~= Censored)+Lam(Y ~= Censored))+ (theta(Y ~= Censored)+Y(Y ~= Censored))./(theta(Y ~= Censored)+Lam(Y ~= Censored));
            gnb(Y ~= Censored,:)  = txnb(Y ~= Censored,:).*gnb_tmp(:, ones(1, NVarNB));
                
            if Truncated == 2
                % Endogenous stratification (on-site sample): the censored cell mixes
                % the size-biased masses P*(k) = k*f(k)/lambda, not f(k). Both the
                % weights and the per-k derivative w.r.t. beta change accordingly:
                %   d log P*(k) / d beta = u*(k-lambda) - 1
                %   d log P*(k) / d s    = d log f(k) / d s   (log k - log lambda is free of gamma)
                pk = fk.*vec./clam;

                % Y<C
                g(Y ~= Censored,:) = g(Y ~= Censored,:) - Xa(Y ~= Censored,:);

                % Y=C: d log ftmp / d beta = [ sum_k P*(k) + u*sum_k P*(k)*(lambda-k) ] / ftmp
                gtmp = (sum(pk,2) + u(:,1).*sum(pk.*(clam-vec),2))./ftmp;
                g(Y == Censored,:) = gtmp(:, ones(1, EstimOpt.NVarA)).*Xa(Y == Censored,:);

                % Y=C: d log ftmp / d s = theta * sum_k P*(k)*dlnf_ds / ftmp
                gtmp = ctheta(:,1).*sum(pk.*dlnf_ds,2)./ftmp;
                gnb(Y == Censored,:) = gtmp(:, ones(1, EstimOpt.NVarNB)).*Xnb(Y == Censored,:);

            else % Truncated == 0 or 1 - the censored cell mixes the plain NB masses
                % d log ftmp / d beta = u * sum_k f(k)*(lambda-k) / ftmp
                gtmp = u(:,1).*sum(fk.*(clam-vec),2)./ftmp;
                g(Y == Censored,:) = gtmp(:, ones(1, EstimOpt.NVarA)).*Xa(Y == Censored,:);

                % d log ftmp / d s = theta * sum_k f(k)*dlnf_ds / ftmp
                gtmp = ctheta(:,1).*sum(fk.*dlnf_ds,2)./ftmp;
                gnb(Y == Censored,:) = gtmp(:, ones(1, EstimOpt.NVarNB)).*Xnb(Y == Censored,:);

                if Truncated == 1
                    % The zero-truncation correction -log(1-u^theta) has the same form
                    % for both groups, so it is subtracted from every row.
                    gtrunc = (ThetaLam.^(theta+1))./(-expm1(theta.*log(ThetaLam)));
                    gtrunc = gtrunc.*Xa.*Lam;
                    g = g - gtrunc;

                    gtruncnb = (ThetaLam.^theta).*(Lam + (theta+Lam).*log(ThetaLam));
                    gtruncnb = gtruncnb./((theta+Lam).*(-expm1(theta.*log(ThetaLam))));
                    gtruncnb = gtruncnb.*Xnb.*theta;
                    gnb = gnb - gtruncnb;
                end
            end
        end


    g = -[g, gnb];
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



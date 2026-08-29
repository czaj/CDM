function nFail = test_CDM_gradients(verbose)
% test_CDM_gradients compares the analytical gradient of LL_CDM against central
% differences over the full grid of NB x NBtype x Zinf x Truncated x Censored x
% RealMin, where NBtype is 2 (NB2), 1 (NB1) or 3 (NB-P).
%
% Returns the number of unexpected failures. Combinations listed in KNOWN below
% are reported but not counted - they document open issues in the audit.

if nargin < 1 || isempty(verbose)
    verbose = 1;
end
cdm_test_paths;

% Open issues - reported as KNOWN, do not fail the suite.
% {NB, NBtype, Zinf, Truncated, Censored, id, note}
% Empty: every combination on the grid matches the numerical gradient, including
% NB1 and NB-P with censoring.
KNOWN = cell(0,7);

n = 2000;
NVarA = 3;
tol = 1e-4;

nFail = 0;
nKnown = 0;
nPass = 0;

if verbose
    fprintf('%-4s %-5s %-5s %-4s %-4s %-5s %13s %12s   %s\n', ...
        'NB','type','Zinf','Trun','Cens','RMin','max rel. err','sum LL','status');
end

for NBpair = [0 2; 1 2; 1 1; 1 3]' % NB and NBtype; NBtype is irrelevant for Poisson
    NB = NBpair(1); NBtype = NBpair(2);
    for Zinf = [0 1]
        for Truncated = [0 1 2]
            for Censored = [0 6]
                for RealMin = [0 1]

                    % zero inflation with zero truncation is not identified - skip
                    if Zinf == 1 && Truncated > 0
                        continue
                    end

                    rng(20260827); % same data for every cell of the grid

                    Xa = [ones(n,1),unifrnd(-1,1,n,1),unifrnd(0,2,n,1)];
                    beta = [1;-0.7;0.4];
                    Lam = exp(Xa*beta);

                    EstimOpt = struct('NB',NB,'Zinf',Zinf,'Truncated',Truncated, ...
                        'Censored',Censored,'NVarA',NVarA,'RealMin',RealMin,'NBtype',NBtype);

                    b = beta;
                    if NB == 1
                        Xnb = ones(n,1);
                        b = [b;0.2];
                        if NBtype == 3
                            b = [b;1.5]; % the power P is estimated
                        end
                        if NBtype == 1
                            th = Lam/0.8; % NB1: theta = lambda/alpha
                        elseif NBtype == 3
                            Ptrue = 1.5;
                            th = Lam.^(2-Ptrue)/0.8; % NB-P: theta = lambda^(2-P)/alpha
                        else
                            th = 2*ones(n,1); % NB2: theta constant
                        end
                        Y = nbinrnd(th,th./(th+Lam));
                    else
                        Xnb = zeros(n,0);
                        Y = poissrnd(Lam);
                    end
                    EstimOpt.NVarNB = size(Xnb,2);

                    if Zinf == 1
                        Xzinf = [ones(n,1),unifrnd(-1,1,n,1)];
                        b = [b;0.3;-0.2];
                        Y = Y.*binornd(1,0.7,n,1);
                    else
                        Xzinf = zeros(n,0);
                    end
                    EstimOpt.NVarZinf = size(Xzinf,2);

                    if Truncated > 0
                        Y = max(Y,1);
                    end
                    if Censored > 0
                        Y = min(Y,Censored);
                    end

                    [isKnown,auditId,note] = lookupKnown(KNOWN,NB,NBtype,Zinf,Truncated,Censored);
                    status = 'ok';
                    relErr = NaN;
                    sumLL = NaN;

                    try
                        [f,g] = LL_CDM(Y,Xa,Xnb,Xzinf,EstimOpt,b);
                        sumLL = -sum(f);
                        ga = sum(g,1)';

                        gn = zeros(numel(b),1);
                        for k = 1:numel(b)
                            h = 1e-6*max(1,abs(b(k)));
                            bp = b; bp(k) = bp(k)+h;
                            bm = b; bm(k) = bm(k)-h;
                            gn(k) = (sum(LL_CDM(Y,Xa,Xnb,Xzinf,EstimOpt,bp)) ...
                                   - sum(LL_CDM(Y,Xa,Xnb,Xzinf,EstimOpt,bm)))/(2*h);
                        end
                        relErr = max(abs(ga-gn)./max(1,abs(gn)));

                        if ~all(isfinite(f))
                            status = 'LL NOT FINITE';
                        elseif relErr > tol
                            status = 'GRAD MISMATCH';
                        end
                    catch ME
                        status = ['ERROR: ' ME.message];
                    end

                    if strcmp(status,'ok')
                        nPass = nPass + 1;
                        if isKnown
                            status = ['ok (' auditId ' appears fixed - update KNOWN)'];
                        end
                    elseif isKnown
                        nKnown = nKnown + 1;
                        status = ['KNOWN ' auditId ': ' status];
                    else
                        nFail = nFail + 1;
                        status = ['FAIL: ' status];
                    end

                    if verbose
                        fprintf('%-4d %-5d %-5d %-4d %-4d %-5d %13.3e %12.4f   %s\n', ...
                            NB,NBtype,Zinf,Truncated,Censored,RealMin,relErr,sumLL,status);
                        if isKnown && contains(status,'KNOWN')
                            fprintf('%56s %s\n','',note);
                        end
                    end

                end
            end
        end
    end
end

if verbose
    fprintf('\ntest_CDM_gradients: %d passed, %d known issues, %d failed\n',nPass,nKnown,nFail);
end
if nFail == 0 && verbose
    disp('test_CDM_gradients passed')
end
end


function [isKnown,auditId,note] = lookupKnown(KNOWN,NB,NBtype,Zinf,Truncated,Censored)
isKnown = false;
auditId = '';
note = '';
for i = 1:size(KNOWN,1)
    if KNOWN{i,1} == NB && KNOWN{i,2} == NBtype && KNOWN{i,3} == Zinf && ...
            KNOWN{i,4} == Truncated && KNOWN{i,5} == Censored
        isKnown = true;
        auditId = KNOWN{i,6};
        note = KNOWN{i,7};
        return
    end
end
end

function nFail = test_CDM_gradients(verbose)
% test_CDM_gradients compares the analytical gradient of LL_CDM against central
% differences over the full grid of NB x Zinf x Truncated x Censored x RealMin.
%
% Returns the number of unexpected failures. Combinations listed in KNOWN below
% are reported but not counted - they document open issues in the audit.

if nargin < 1 || isempty(verbose)
    verbose = 1;
end
cdm_test_paths;

% Open issues - reported as KNOWN, do not fail the suite.
% {NB, Zinf, Truncated, Censored, audit id, note}
% Empty: every combination on the grid currently matches the numerical gradient.
KNOWN = cell(0,6);

n = 2000;
NVarA = 3;
tol = 1e-4;

nFail = 0;
nKnown = 0;
nPass = 0;

if verbose
    fprintf('%-6s %-5s %-4s %-4s %-5s %13s %12s   %s\n', ...
        'NB','Zinf','Trun','Cens','RMin','max rel. err','sum LL','status');
end

for NB = [0 1]
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
                        'Censored',Censored,'NVarA',NVarA,'RealMin',RealMin);

                    b = beta;
                    if NB == 1
                        Xnb = ones(n,1);
                        b = [b;0.2];
                        Y = nbinrnd(2,2./(2+Lam));
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

                    [isKnown,auditId,note] = lookupKnown(KNOWN,NB,Zinf,Truncated,Censored);
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
                        fprintf('%-6d %-5d %-4d %-4d %-5d %13.3e %12.4f   %s\n', ...
                            NB,Zinf,Truncated,Censored,RealMin,relErr,sumLL,status);
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


function [isKnown,auditId,note] = lookupKnown(KNOWN,NB,Zinf,Truncated,Censored)
isKnown = false;
auditId = '';
note = '';
for i = 1:size(KNOWN,1)
    if KNOWN{i,1} == NB && KNOWN{i,2} == Zinf && KNOWN{i,3} == Truncated && KNOWN{i,4} == Censored
        isKnown = true;
        auditId = KNOWN{i,5};
        note = KNOWN{i,6};
        return
    end
end
end

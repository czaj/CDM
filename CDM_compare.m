function Test = CDM_compare(Results1,Results2,Display)
% CDM_compare tests one estimated count data model against another.
%
% Syntax:   Test = CDM_compare(ResultsRestricted,ResultsUnrestricted)
%           Test = CDM_compare(Results1,Results2,Display)
%
% Which test is run depends on how the two models relate:
%
%   Poisson inside Negative Binomial - a likelihood ratio test on the boundary of
%     the parameter space. Poisson is NB with alpha = 0, which sits on the edge of
%     the admissible set, so the statistic is distributed as an equal mixture of a
%     point mass at zero and chi2(1). Using plain chi2(1) makes the test roughly
%     twice as conservative as it should be.
%
%   zero-inflated against its plain counterpart - these are not nested in the
%     usual sense (the inflation probability is unidentified when there is no
%     excess of zeros), so the Voung test for non-nested models is used instead.
%     It needs Results.LLdetailed, which CDM always stores.
%
%   anything else - a likelihood ratio test with chi2(dk), reported with a warning
%     that nesting is assumed rather than checked.
%
% See also CDM, CDM_post
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

if nargin < 3 || isempty(Display)
    Display = 1;
end

E1 = Results1.EstimOpt;
E2 = Results2.EstimOpt;
LL1 = Results1.LL;
LL2 = Results2.LL;
k1 = E1.Params;
k2 = E2.Params;

Test = struct('type','','statistic',NaN,'df',NaN,'pValue',NaN,'conclusion','');

sameCount = isequal(E1.Truncated,E2.Truncated) && isequal(E1.Censored,E2.Censored);

if sameCount && E1.Zinf == E2.Zinf && E1.NB == 0 && E2.NB == 1
    % Poisson nested in NB, on the boundary alpha = 0
    LR = 2*(LL2-LL1);
    LR = max(LR,0);
    p = 0.5*(1-chi2cdf(LR,1)); % 0.5*chi2(0) + 0.5*chi2(1)
    Test.type = 'LR, Poisson vs Negative Binomial (boundary-corrected)';
    Test.statistic = LR; Test.df = 1; Test.pValue = p;
    if p < 0.05
        Test.conclusion = 'overdispersion: reject Poisson in favour of NB';
    else
        Test.conclusion = 'no evidence of overdispersion: Poisson is adequate';
    end

elseif sameCount && E1.NB == E2.NB && E1.Zinf ~= E2.Zinf
    % zero inflation is not identified under the null - Voung rather than LR
    if ~isfield(Results1,'LLdetailed') || ~isfield(Results2,'LLdetailed')
        error('CDM_compare: the Voung test needs Results.LLdetailed from both models')
    end
    if E1.Zinf == 1 % the plain model has to go first for the sign convention below
        [Results1,Results2] = deal(Results2,Results1);
        [k1,k2] = deal(k2,k1);
    end
    Test.type = 'Voung, plain vs zero-inflated (non-nested)';
    V = voung_test(Results1.LLdetailed,Results2.LLdetailed,k1,k2); % returns [statistic, p]
    Test.statistic = V(1);
    Test.pValue = V(2);
    if Test.statistic < -1.96
        Test.conclusion = 'favours the zero-inflated model';
    elseif Test.statistic > 1.96
        Test.conclusion = 'favours the plain model';
    else
        Test.conclusion = 'neither model is favoured';
    end

else
    dk = k2-k1;
    if dk <= 0
        error('CDM_compare: pass the restricted model first - it must have fewer parameters')
    end
    LR = max(2*(LL2-LL1),0);
    Test.type = 'LR (nesting assumed, not verified)';
    Test.statistic = LR; Test.df = dk;
    Test.pValue = 1-chi2cdf(LR,dk);
    Test.conclusion = 'check that the models really are nested';
end

if Display ~= 0
    fprintf('\n')
    cprintf('*Black',['Model comparison: ',Test.type,' \n'])
    fprintf('%-26s%12.4f\n','LL restricted:',min(LL1,LL2))
    fprintf('%-26s%12.4f\n','LL unrestricted:',max(LL1,LL2))
    fprintf('%-26s%12.4f\n','statistic:',Test.statistic)
    if ~isnan(Test.df)
        fprintf('%-26s%12d\n','df:',Test.df)
    end
    fprintf('%-26s%12.4f\n','p-value:',Test.pValue)
    fprintf('%-26s%s\n','conclusion:',Test.conclusion)
    fprintf('\n')
end
end

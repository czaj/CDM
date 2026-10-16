function [LL,grad,hess] = LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOpt,b)
% LL_CDM_hybrid_MATlike Negative joint likelihood for MATLAB optimizers.
%   [LL,grad] = LL_CDM_hybrid_MATlike(INPUT,EstimOpt,OptimOpt,b)
%   uses resolved inputs/options from DataCleanCDM_hybrid. OptimOpt is accepted
%   for the standard CDM signature; the likelihood does not depend on it.
%   The gradient is analytic. Optional hess is the BHHH score cross-product,
%   NOT the observed Hessian used by CDM_hybrid for inference. Supplied weights
%   enter the likelihood and scores once; BHHH consequently squares weights.
%
%   See also LL_CDM_hybrid, CDM_hybrid, LL_CDM_MATlike.
[lli,score] = LL_CDM_hybrid(b,INPUT,EstimOpt);
LL = -sum(lli);
grad = -sum(score,1)';
if nargout>2, hess = score'*score; end
if ~isfinite(LL) || any(~isfinite(grad))
    % Reject an unrepresentable optimizer trial without clipping parameters
    % or replacing small probabilities by a likelihood-changing floor.
    LL = Inf;
    grad = zeros(size(b));
    if nargout>2, hess = nan(numel(b)); end
end
end

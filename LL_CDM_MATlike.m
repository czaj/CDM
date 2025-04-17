% function [LL,g,h]= LL_CDM_MATlike(Y,Xa,Xnb,W,EstimOpt,OptimOpt,b0)
function [LL,g,h]= LL_CDM_MATlike(Y,Xa,Xnb,Xzinf,W,EstimOpt,OptimOpt,b0)

% LLfun = @(B) LL_CDM(Y,Xa,Xnb, EstimOpt,B);
%powinno być
LLfun = @(B) LL_CDM(Y,Xa,Xnb, Xzinf, EstimOpt,B);


if isequal(OptimOpt.GradObj,'on')
    if EstimOpt.NumGrad == 0
        [f,j] = LLfun(b0);
        j(:,EstimOpt.BActive ==0) = 0;
        f = f.*W;
        j = j.*W(:, ones(1, size(j,2)));
        g = sum(j,1)'; ...
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    else % elseif EstimOpt.NumGrad == 1 
        f = LLfun(b0);  
        j = numdiff(LLfun,f,b0,isequal(OptimOpt.FinDiffType,'central'),EstimOpt.BActive);...
        f = f.*W;
        j = j.*W(:, ones(1, size(j,2)));
        g = sum(j,1)';   
        if isequal(OptimOpt.Hessian,'user-supplied') == 1
            h = j'*j;
        end
    end
    
else % No gradient
    EstimOpt.NumGrad = 1;
    f = LLfun(b0).*W;
end
LL = sum(f);

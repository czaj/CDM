function Report = CDM_checkdesign(INPUT,EstimOpt)
% CDM_checkdesign screens the explanatory variables before estimation.
%
% Syntax:   Report = CDM_checkdesign(INPUT,EstimOpt)
%
% Called by CDM and CDM_hurdle. There is normally no reason to call it directly,
% but it is a separate function so that it can be tested and so that a user who
% wants the diagnostics without estimating anything can have them.
%
% WHY THIS EXISTS
%   A count model whose design matrix is rank deficient still converges. The
%   log-likelihood is the log-likelihood of the same fitted model, the gradient
%   is zero at the optimum, and the reported fit statistics are correct. What is
%   not correct is any individual coefficient: the optimiser splits the shared
%   direction arbitrarily between the collinear columns, and where it lands
%   depends on the starting values.
%
%   For most models that is a nuisance. For a travel cost model it is fatal,
%   because the headline output is consumer surplus, -1/beta, read from one
%   coefficient slot. Duplicating the travel cost column in an otherwise sound
%   design leaves the log-likelihood identical to three decimals and changes the
%   reported surplus from 1.567 to 47.978. Near-duplication is worse: it does
%   not even produce the NaN standard error that exact duplication does, and at
%   r = 0.999 the same surplus comes out as 0.035. Nothing downstream notices,
%   because nothing downstream has any way to.
%
% WHAT IT DOES
%   rank            An exactly rank-deficient design is an error. The parameters
%                   are not identified and no amount of care downstream can
%                   recover them. The columns responsible are named, found by a
%                   QR with column pivoting. Set EstimOpt.AllowRankDeficient = 1
%                   to downgrade this to a warning, which is occasionally what
%                   you want when parameters are being constrained with BActive.
%
%   conditioning    The condition index of the column-scaled design. Above 30 is
%                   the conventional threshold for collinearity worth worrying
%                   about (Belsley, Kuh and Welsch); above 100 the coefficients
%                   should not be interpreted individually. A warning, not an
%                   error, because a badly conditioned design is still estimable
%                   and the user may know exactly what they are doing.
%
%   travel cost     The variance inflation factor of the parameter named by
%                   EstimOpt.TC_par_index, which is the one consumer surplus is
%                   read from. A VIF above 10 is reported with the consequence
%                   spelled out, because that coefficient is not just another
%                   coefficient here.
%
%   Xnb, Xzinf      The same rank and conditioning checks, since a collinear
%                   overdispersion or inflation equation breaks in the same way.
%
% RETURNS  a struct with one field per design matrix, each holding the rank, the
% number of columns, the condition index and the names of any dependent columns,
% plus tcVIF for the travel cost parameter. Nothing is printed when
% EstimOpt.Display is 0.
%
% See also CDM, CDM_hurdle, DataCleanCDM
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

if nargin < 2
    error('Too few input arguments for CDM_checkdesign(INPUT,EstimOpt)')
end
if ~isfield(EstimOpt,'Display') || isempty(EstimOpt.Display)
    EstimOpt.Display = 1;
end
if ~isfield(EstimOpt,'AllowRankDeficient') || isempty(EstimOpt.AllowRankDeficient)
    EstimOpt.AllowRankDeficient = 0;
end
Display = EstimOpt.Display;

Report = struct();

blocks = {'Xa','NamesA','count equation'; ...
          'Xnb','NamesNB','overdispersion equation'; ...
          'Xzinf','NamesZinf','zero-inflation equation'};

for b = 1:size(blocks,1)
    fld = blocks{b,1};
    if ~isfield(INPUT,fld) || isempty(INPUT.(fld)) || size(INPUT.(fld),2) == 0
        continue
    end
    X = INPUT.(fld);
    names = local_names(EstimOpt, blocks{b,2}, size(X,2));

    if any(~isfinite(X(:)))
        error('INPUT.%s contains non-finite values (%d of %d elements)', ...
            fld, sum(~isfinite(X(:))), numel(X))
    end

    R = local_one(X, names, fld, blocks{b,3}, Display, EstimOpt.AllowRankDeficient);
    Report.(fld) = R;
end

% The travel cost parameter is not just another coefficient: consumer surplus is
% -1/beta read from that slot, so a coefficient that is only weakly separated
% from the others turns into a surplus that is mostly noise.
Report.tcVIF = NaN;
if isfield(EstimOpt,'TC_par_index') && ~isempty(EstimOpt.TC_par_index) && ...
        isfield(INPUT,'Xa') && ~isempty(INPUT.Xa)
    j = EstimOpt.TC_par_index;
    if j >= 1 && j <= size(INPUT.Xa,2)
        Report.tcVIF = local_vif(INPUT.Xa, j);
        if Display ~= 0 && isfinite(Report.tcVIF) && Report.tcVIF > 10
            nm = local_names(EstimOpt,'NamesA',size(INPUT.Xa,2));
            local_warn(sprintf(['the travel cost parameter (%s) has a variance ' ...
                'inflation factor of %.1f.\n         Consumer surplus is -1/beta ' ...
                'read from that one coefficient, so most of what\n         is ' ...
                'reported as a surplus here is the arbitrary split of a shared ' ...
                'direction.\n'], nm{j}, Report.tcVIF));
        end
    end
end
end


% ======================================================================= utils
function R = local_one(X, names, fld, what, Display, allowDeficient)
k = size(X,2);
R = struct('rank',k,'ncol',k,'conditionIndex',NaN,'dependent',{{}});

% --- rank, with the offending columns named.  QR with column pivoting orders
% the columns by how much they add, so anything past the numerical rank is
% linearly dependent on what came before it.
[~, Rq, E] = qr(X, 0);
dR = abs(diag(Rq));
tol = max(size(X)) * eps(max(dR));
r = sum(dR > tol);
R.rank = r;
if r < k
    R.dependent = names(sort(E(r+1:end)));
    msg = sprintf(['INPUT.%s is rank deficient: rank %d of %d columns.\n' ...
        '         Linearly dependent on the others: %s\n' ...
        '         The likelihood will still converge and the fit statistics will\n' ...
        '         still be right, but the individual coefficients are not\n' ...
        '         identified - the optimiser splits the shared direction between\n' ...
        '         the collinear columns arbitrarily.  Consumer surplus is read\n' ...
        '         from one such coefficient, so it would be reported as a number\n' ...
        '         with no meaning.  Drop the redundant column(s), or set\n' ...
        '         EstimOpt.AllowRankDeficient = 1 if the design is deliberate.'], ...
        fld, r, k, strjoin(R.dependent, ', '));
    if allowDeficient
        local_warn([msg newline]);
    else
        error('CDM:RankDeficientDesign', '%s', msg);
    end
end

% --- conditioning of the column-scaled design
nrm = sqrt(sum(X.^2,1));
nrm(nrm == 0) = 1;
Xs = X ./ nrm;
s = svd(Xs);
if ~isempty(s) && s(end) > 0
    R.conditionIndex = s(1)/s(end);
end
if Display ~= 0 && isfinite(R.conditionIndex) && R.conditionIndex > 30 && r == k
    if R.conditionIndex > 100
        local_warn(sprintf(['the %s has a condition index of %.0f.\n' ...
            '         Above 100 the coefficients should not be read individually;\n' ...
            '         the model as a whole is still fine.\n'], what, R.conditionIndex));
    else
        local_warn(sprintf(['the %s has a condition index of %.0f (above 30).\n' ...
            '         Collinearity is inflating the standard errors.\n'], ...
            what, R.conditionIndex));
    end
end
end

function v = local_vif(X, j)
% Variance inflation factor for column j: regress it on the rest of the design.
% A constant already in X does the work of an intercept.
k = size(X,2);
if k < 2; v = 1; return; end
y = X(:,j);
Z = X(:,[1:j-1, j+1:k]);
if all(abs(y - y(1)) < eps(max(abs(y),1))); v = NaN; return; end   % constant column
if rcond(Z'*Z) < eps; v = Inf; return; end
res = y - Z*(Z\y);
sst = sum((y - mean(y)).^2);
if sst <= 0; v = NaN; return; end
r2 = 1 - sum(res.^2)/sst;
r2 = min(max(r2,0), 1 - eps);
v = 1/(1 - r2);
end

function names = local_names(EstimOpt, field, k)
names = cell(k,1);
if isfield(EstimOpt, field) && ~isempty(EstimOpt.(field)) && numel(EstimOpt.(field)) == k
    src = EstimOpt.(field);
    for i = 1:k
        if iscell(src); names{i} = char(string(src{i})); else; names{i} = char(string(src(i))); end
    end
else
    for i = 1:k; names{i} = sprintf('column %d', i); end
end
end

function local_warn(msg)
if exist('cprintf','file') == 2 && exist('rgb','file') == 2
    cprintf(rgb('DarkOrange'), ['WARNING: ' msg]);
else
    fprintf(2, ['WARNING: ' msg]);
end
end

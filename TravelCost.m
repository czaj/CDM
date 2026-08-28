function [TC,Detail] = TravelCost(Trip,Opt)
% TravelCost builds the travel cost variable for a travel cost study.
%
% Syntax:   [TC,Detail] = TravelCost(Trip)
%           [TC,Detail] = TravelCost(Trip,Opt)
%
% Travel cost is the money spent getting to the site plus the value of the time
% spent doing so. The money part is uncontroversial; the time part is not, which
% is why the rule is an explicit choice here and why the default returns a whole
% sensitivity panel rather than one number.
%
% Trip:
%    Distance     distance to the site, in km. One-way unless Trip.RoundTrip = 1
%    CostPerKm    monetary cost per km (scalar, or one value per respondent)
%    Time         travel time, in hours, on the same one-way/round-trip basis
%    Wage         hourly wage or income rate (scalar, or one value per respondent)
%    OnSiteTime   optional time spent on site, in hours, valued by the same rule
%    ValuesTime   optional 0/1 indicator, required by the 'full_or_zero' rule:
%                 1 for respondents who would have preferred a shorter trip,
%                 0 for those who would not
%    RoundTrip    1 if Distance and Time are already round-trip (default 0)
%
% Opt:
%    TimeRule     how to value travel time:
%       'fraction'      a fraction of the wage - the usual convention. Returns one
%                       column per fraction in Opt.TimeFraction (default
%                       [0 1/4 1/3 1/2]), which is the sensitivity analysis
%                       reviewers ask for
%       'full'          the full wage rate
%       'full_or_zero'  the full wage for respondents who would have preferred a
%                       shorter trip, and nothing for those who would not. Travel
%                       that someone would not shorten is not time given up, so
%                       charging them for it overstates the price they faced
%       'none'          money cost only
%    TimeFraction fraction(s) of the wage under the 'fraction' rule
%    Display      0 to compute quietly
%
% Detail returns the money and time components separately, the rule used and the
% column labels, so the panel can go straight into a table.
%
% Example:
%    Trip = struct('Distance',dist,'CostPerKm',0.35,'Time',dist/70,'Wage',wage);
%    [TC,D] = TravelCost(Trip);            % sensitivity over the usual fractions
%    INPUT.Xa = [ones(n,1), TC(:,3), inc]; % 1/3 of the wage
%
% See also CDM, CDM_post
%
% Author: Mikolaj Czajkowski, Professor
% University of Warsaw, Faculty of Economic Sciences
% email address: mik@czaj.org
% Website: http://czaj.org/#

if nargin < 2 || isempty(Opt)
    Opt = struct();
end
if ~isfield(Opt,'TimeRule') || isempty(Opt.TimeRule)
    Opt.TimeRule = 'fraction';
end
if ~isfield(Opt,'TimeFraction') || isempty(Opt.TimeFraction)
    Opt.TimeFraction = [0 1/4 1/3 1/2];
end
if ~isfield(Opt,'Display')
    Opt.Display = 1;
end

req = {'Distance','CostPerKm','Time','Wage'};
for i = 1:numel(req)
    if ~isfield(Trip,req{i}) || isempty(Trip.(req{i}))
        error('TravelCost: Trip.%s is required',req{i})
    end
end

dist = Trip.Distance(:);
time = Trip.Time(:);
n = numel(dist);
if numel(time) ~= n
    error('TravelCost: Trip.Distance and Trip.Time must have the same length')
end

if ~isfield(Trip,'RoundTrip') || isempty(Trip.RoundTrip) || Trip.RoundTrip == 0
    dist = 2*dist; % there and back
    time = 2*time;
end

costPerKm = expandTo(Trip.CostPerKm,n,'CostPerKm');
wage = expandTo(Trip.Wage,n,'Wage');

if isfield(Trip,'OnSiteTime') && ~isempty(Trip.OnSiteTime)
    timeValued = time + expandTo(Trip.OnSiteTime,n,'OnSiteTime');
else
    timeValued = time;
end

if any(dist < 0) || any(timeValued < 0) || any(wage < 0)
    error('TravelCost: distance, time and wage must be non-negative')
end

Money = dist.*costPerKm;

switch lower(Opt.TimeRule)
    case 'none'
        rate = zeros(n,1);
        labels = {'money only'};
    case 'full'
        rate = wage;
        labels = {'full wage'};
    case 'full_or_zero'
        if ~isfield(Trip,'ValuesTime') || isempty(Trip.ValuesTime)
            error(['TravelCost: the ''full_or_zero'' rule needs Trip.ValuesTime - ' ...
                '1 for respondents who would have preferred a shorter trip, 0 otherwise'])
        end
        v = Trip.ValuesTime(:);
        if numel(v) ~= n
            error('TravelCost: Trip.ValuesTime must have one entry per respondent')
        end
        if ~all(ismember(v,[0 1]))
            error('TravelCost: Trip.ValuesTime must contain only 0 and 1')
        end
        rate = wage.*v;
        labels = {'full wage or zero'};
    case 'fraction'
        fr = Opt.TimeFraction(:)';
        rate = wage*fr; % n x numel(fr)
        labels = arrayfun(@(f) sprintf('%.3g of wage',f),fr,'UniformOutput',false);
    otherwise
        error('TravelCost: unknown Opt.TimeRule ''%s''',Opt.TimeRule)
end

Time = timeValued.*rate;
TC = Money + Time;

Detail = struct('Money',Money,'Time',Time,'Rule',Opt.TimeRule, ...
    'Labels',{labels},'Distance',dist,'TimeHours',timeValued);
if strcmpi(Opt.TimeRule,'fraction')
    Detail.TimeFraction = Opt.TimeFraction(:)';
end
if strcmpi(Opt.TimeRule,'full_or_zero')
    Detail.ShareValuingTime = mean(Trip.ValuesTime(:));
end

if Opt.Display ~= 0
    fprintf('\n')
    cprintf('*Black','Travel cost \n')
    fprintf('%-28s%12.4f\n','money cost, mean:',mean(Money))
    fprintf('%-28s%12.4f\n','round-trip km, mean:',mean(dist))
    fprintf('%-28s%12.4f\n','hours valued, mean:',mean(timeValued))
    if strcmpi(Opt.TimeRule,'full_or_zero')
        fprintf('%-28s%12.4f\n','share whose time is valued:',Detail.ShareValuingTime)
    end
    fprintf('%-28s%12s%12s%12s\n','time valuation','time cost','total','share time')
    for c = 1:size(TC,2)
        fprintf('%-28s%12.4f%12.4f%12.3f\n',labels{c}, ...
            mean(Time(:,c)),mean(TC(:,c)),mean(Time(:,c))/mean(TC(:,c)))
    end
    fprintf('\n')
end
end


function v = expandTo(v,n,name)
v = v(:);
if isscalar(v)
    v = v*ones(n,1);
elseif numel(v) ~= n
    error('TravelCost: Trip.%s must be a scalar or have one value per respondent',name)
end
end

function varargout = cohortplot(data, rank_criteria, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% cohortaxes = cohortplot(data, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% COHORTPLOT takes single-individual trajectories (rows of data_in), ranks them along specified 
% criteria (e.g. summed activity - provided as a vector), and outputs a line plot where trajectories
% are averaged by a percentile "cohort" - similar traces, as determined by ranking criteria. 
% 
% INPUTS (required)
% data                matrix of (e.g. single-cell) trajectories -> size = [n m], where n is the # of individuals
%                     and m is the number of observations (e.g. timepoints)
% ranking_criteria    vector (of size [n x 1]) that are ranked to group data trajectories by similarity
%
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Colors'         [N x 3] matrix specifying line colors - cycles if N < number of cohorts
% 'NumCohorts'     Number of cohorts - e.g. 4 would group 0th-25th, 25th-50th, 50th-75th, and 75th-100th
%                  percentile data. Default is 8.
% 'SortDirection'  Direction to sort ranking - either 'ascend' or 'descend'. NaN values will be omitted.
% 'XVector'        Alternate x (time) vector to use
% 'LineWidth'      Plot line width
% 'Legend'         Specifies legend location (or 'none', if no legend is desired)
% 'Axes'           Axes handle where graph will be plotted (default: create new figure + axes)
%
%
% OUTPUTS
% cohort_axes        Axes handle of cohort plot figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% 1) Plot data (matrix)
valid_data = @(x) assert(ismatrix(x) && isnumeric(x),'Plot data must be numerical matrix');
addRequired(p,'data',valid_data);
% 2) Rank data (vector)
valid_rank = @(x) assert(ismatrix(x) && numel(x)==size(data,1),'Plot data must be numerical matrix');
addRequired(p,'rank_criteria',valid_rank);

% Optional parameters
valid_color = @(x) assert(isempty(x)||(ismatrix(x)&&size(x,2)==3), 'Specify colors with a matrix of RGB triplets (similar to colormap)');
addParameter(p,'Color', [],valid_color);
addParameter(p,'NumCohorts',8,@isnumeric)
valid_dir = @(x) assert(strcmp(x,'ascend')||strcmp(x,'descend'),'Sort order must be either ''ascend'' or ''descend''');
addParameter(p,'SortDirection','ascend',valid_dir) 
valid_vect = @(x) assert(isnumeric(x) && numel(x)==size(data,2),'X vector must be same size as # of cols in input data');
addParameter(p,'XVector',1:size(data,2),valid_vect) 
addParameter(p,'LineWidth', 1, @isnumeric);
valid_legend = @(x) assert(ischar(x),'Specify either a valid location (e.g. ''northeast'') or ''none''');
addParameter(p,'Legend', 'none', valid_legend);
addParameter(p,'Axes',nan,@ishandle);

% Parse inputs, save some to variables
parse(p,data, rank_criteria, varargin{:})
clr = p.Results.Color;
cohorts = p.Results.NumCohorts;
t = p.Results.XVector;


%% Set line colors
if isempty(clr)
    colormaps = loadcolormaps;
    clr = colormaps.magma(round(linspace(1,220,cohorts)),:);
end
if size(clr,1)<cohorts
    clr = repmat(clr,[ceil(cohorts/size(clr,1)) 1]);
end

%% Rank data, do averaging across cohorts
if numel(rank_criteria) ~= size(data,1)
    error('Number of elements in ranking criteria, rank_criteria, must correspond to number of rows in data')
end

data(isnan(rank_criteria),:) = [];
rank_criteria(isnan(rank_criteria),:) = [];
[rank_criteria, order] = sort(rank_criteria,'ascend');
data = data(order,:);
groups = floor(linspace(1,length(data),cohorts+1));
pctiles = floor(linspace(1,100,cohorts+1));

cohort_trajectories = zeros(cohorts,size(data,2));
legend_entries = cell(1,cohorts);
for i = 1:cohorts
    cohort_trajectories(i,:) = nanmean(data(groups(i):groups(i+1),:));
    if i==1
        legend_entries{i} = ['1st - ',num2ordinal(pctiles(i+1),0), ' pctile'];
    else
        legend_entries{i} = [num2ordinal(pctiles(i)+1,0),' - ',num2ordinal(pctiles(i+1),0), ' pctile'];

    end
end




%% 

% Create figure (if axes wasn't provided)
if ~ishandle(p.Results.Axes)
    cohort_axes = axes('Parent',gcf);
else
    cohort_axes = p.Results.Axes;
end

hold(cohort_axes,'on')
set(cohort_axes,'ColorOrder',clr)

plot(t(:), cohort_trajectories', 'LineWidth', p.Results.LineWidth)
hold(cohort_axes,'off')

if ~strcmp(p.Results.Legend,'none')
    legend(legend_entries,'Location',p.Results.Legend)
end

if nargout > 0
    varargout{1} = cohort_axes;
end


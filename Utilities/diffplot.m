function varargout = diffplot(data,rank_criteria, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% diff_fig = diffplot(data, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DIFFPLOT pairs a "cohort plot" (which averages progression of cohorts of single trajectories, grouped 
% into bins by state variable (e.g. end-state expression of a protein). If a threshold is specified,
% cohorts will be colored according to whether they exceed this threshold or not.
% 
% INPUTS (required)
% data                matrix of (e.g. single-cell) trajectories -> size = [n m], where n is the # of individuals
%                     and m is the number of observations (e.g. timepoints)
% ranking_criteria    vector (of size [n x 1]) that are ranked to group data trajectories by similarity
%
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Color'         [N x 3] matrix specifying line colors - cycles if N < number of cohorts
% 'NumCohorts'     Number of cohorts - e.g. 4 would group 0th-25th, 25th-50th, 50th-75th, and 75th-100th
%                  percentile data. Default is 8.
% 'Threshold'      Single value used to divide cells (e.g. differentiated vs undifferentiated)
% 'XVector'        Alternate x (time) vector to use
% 'LineWidth'      Plot line width
% 'Legend'         Specifies legend location (or 'none', if no legend is desired)
% 'Bins'           Specify bins (or bin limits) of histogram.

%
% OUTPUTS
% cohort_handle    Axes handle of cohort plot (passed as varargout{1})
% hist_handle    Axes handle of cohort plot (passed as varargout{2})

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
addParameter(p,'Threshold',nan,@isnumeric)
valid_vect = @(x) assert(isnumeric(x) && numel(x)==size(data,2),'X vector must be same size as # of cols in input data');
addParameter(p,'XVector',1:size(data,2),valid_vect) 
addParameter(p,'LineWidth', 1, @isnumeric);
valid_legend = @(x) assert(ischar(x),'Specify either a valid location (e.g. ''northeast'') or ''none''');
addParameter(p,'Legend', 'none', valid_legend);
addParameter(p,'Axes',nan,@ishandle);
valid_bins = @(x) assert(isnumeric(x)&&length(x)>=2,'Bins vector should be numeric and either specify all bins, or [min max]');
addParameter(p,'Bins',nan, @isnumeric);



% Parse inputs, save some to variables
parse(p,data, rank_criteria, varargin{:})
clr = p.Results.Color;
numcohorts = p.Results.NumCohorts;
t = p.Results.XVector;
bins = p.Results.Bins;
thresh1 = p.Results.Threshold;

if isempty(clr)
    clr = [0.9059 0.2980 0.2353];
end

% Filter out any trajectories w/ invalid ranking criteria
data(isnan(rank_criteria),:) = [];
rank_criteria(isnan(rank_criteria),:) = [];


% If a threshold is defined, define the cohort colors accordingly. If < 2 were passed, add gray 
if ~isnan(thresh1)  
    if size(clr,1)<2
        clr = [0.6275 0.6471 0.6667; clr];
    end
    pct_lo = mean(rank_criteria<thresh1);
    clr = [repmat(clr(1,:),[round(numcohorts*pct_lo) 1]); repmat(clr(2,:),[round(numcohorts*(1-pct_lo)) 1]) ];    
end


figure('Position',positionfig(736,350));

% LEFT plot: cohort progression.
h1 = subplot(1,4,1:3);
cohortplot(data,nanmean(data(:,end-20:end-10),2),'NumCohorts',numcohorts,'Axes',h1,'LineWidth',p.Results.LineWidth,...
    'Color',clr,'XVector',p.Results.XVector,'Legend',p.Results.Legend)
set(h1,'XLim',[min(t) max(t)],'FontSize',12);


% RIGHT plot: end-state histogram. Color by threshold, if defined.
h2 = subplot(1,4,4);
if isnan(bins)
    bins = linspace(min(get(h1,'YLim')),max(get(h1,'YLim')),64); % Define bins using plot limits
elseif length(bins)==2
    bins = linspace(min(bins), max(bins), 64); % Define bins 
end

h = histcounts(nanmean(rank_criteria,2),bins);
bin_centers = bins(1:end-1) + diff(bins)/2;
if isnan(thresh1)
barh(bin_centers,h,1,'Parent',h2,'FaceColor',[0.6275 0.6471 0.6667],'EdgeColor','none')
else
    hold(h2,'on')
    barh(bin_centers(bin_centers<thresh1),h(bin_centers<thresh1),1,'Parent',h2,'FaceColor',clr(1,:),'EdgeColor','none')
    barh(bin_centers(bin_centers>=thresh1),h(bin_centers>=thresh1),1,'Parent',h2,'FaceColor',clr(end,:),'EdgeColor','none')
    hold(h2,'off')
end
    
set(h2,'YLim',get(h1,'YLim'),'YTick',get(h1,'YTick'),'YTickLabel',{},'Box','off','FontSize',12)



if nargout > 0
    varargout{1} = h1;
    if nargout > 1
        varargout{2} = h2;
    end
end

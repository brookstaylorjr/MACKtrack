function swarmaxes = beeswarm(data, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% handles = beeswarm(data, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BEESWARM plots distributions of points by spreading them around the y-axis. Input data
% can be 1) a matrix of column data, or (2) a cell matrix (one set of "swarms" per cell - will
% be subgrouped by columns in the cell matrix)
% 
% INPUTS (required)
% data           matrix or cell array of 
%
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Color'        [N x 3] matrix specifying fill colors for points - cycles if length < N
% 'RowPoints     Number of points put in a "row" of the swarm - increase to shrink spacing between points (default = 15)
% 'RowNum'       Number of rows/bins in swarm - increase to place more points along centerline (default = 12)
% 'SwarmWidth'   Relative width of each swarm (scaled from 1.0). Increase to make swarms wider/closer together
% 'GroupSpacing' Spacing between separate groups of swarms (which may be composed of subgroups)
% 'MarkerSize'   Size of each plotted point (default = 7)
% 'ShowMeans'    Show bar for mean. Default = TRUE
% 'Axes'         Axes handle where graph will be plotted (default: create new figure + axes)
%
%
% OUTPUTS
% swarm_axes        Axes handle of violin figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% 1) Plot data (matrix or cell)
valid_data = @(x) assert(ismatrix(x) || iscell(x),'Plot data amust be matrix or cell matrix');
addRequired(p,'data',valid_data);

% Optional parameters
valid_color = @(x) assert(isempty(x)||(ismatrix(x)&&size(x,2)==3), 'Specify colors with a matrix of RGB triplets (similar to colormap)');
addParameter(p,'Color', [],valid_color);
%addParameter(p,'Labels', {},@iscell);
addParameter(p,'RowPoints',15,@isnumeric)
addParameter(p,'RowNum',12,@isnumeric)
addParameter(p,'SwarmWidth',1,@isnumeric)
addParameter(p,'GroupSpacing',1,@isnumeric)
addParameter(p,'MarkerSize',7,@isnumeric)
addParameter(p,'ShowMeans',1,@isboolean)
addParameter(p,'Axes',nan,@ishandle);

% Parse inputs, save some to variables
parse(p,data, varargin{:})
clr = p.Results.Color;
n_pts = p.Results.RowPoints;
n_bins = p.Results.RowNum;
swarm_width = p.Results.SwarmWidth;
marker_size = p.Results.MarkerSize;
show_means = p.Results.ShowMeans;
group_space = p.Results.GroupSpacing;


% RESTRUCTURE data - should be cell matrix
if ~iscell(data)
    tmp = data;
    data = cell(1,1);
    data{1} = tmp;
else
    data = data(:)'; % Flatten cells
    
end

% SET SPACING for groups and subgroups
% Subgroup spacing is on integral values (e.g. 1,2,3...). Spacing between groups is scaled by group size
sz_func = @(x) size(x,2);
n_subgroups = max(cellfun(sz_func,data));
base_spacing = 1:size(data{1},2);
if length(data)>1
    for i = 2:length(data)
        base_spacing = [base_spacing, max(base_spacing) + n_subgroups/2*group_space + (1:size(data{i},2))];
    end
end

num = n_subgroups;

if num==1; num = length(data); end
if isempty(clr)
    clr = cbrewer('div','Spectral', num);
end
% Cycle colors, if required
if size(clr,1)<num
    clr = repmat(clr,ceil(nim/size(clr,1)),1);
end

% SET SPACING for individual points within a subgroup
all_x = data;
if mod(n_pts,2)==0; n_pts = n_pts+1; end % Make sure n_pts is odd
odd_spacings = linspace(-1,1,n_pts);
even_spacings = linspace(-1,1,n_pts+1);

idx = 0;
for i = 1:length(data)
    for j = 1:size(data{i},2)
        col = data{i}(:,j);
        idx= idx+1;
        % Bin points from 1st-99th percentile.
        bin_edges = linspace(prctile(col,1),prctile(col,99),n_bins);
        bin_edges(1) = min(col); bin_edges(end) = max(col);
        [~,~, bin] = histcounts(col,bin_edges);
        x_vals = ones(size(col));
        % Space out any points within the same bin.
        for k = 1:n_bins    
            bin_count = sum(bin==k);
            if bin_count > (n_pts+1)
                spacings = linspace(-1,1,bin_count);
            else
                if mod(bin_count,2) == 0 % Even
                    spacings = even_spacings( (((n_pts+1)/2) - (bin_count/2) + 1) : (((n_pts+1)/2) + (bin_count/2)));
                else % Odd
                    spacings = odd_spacings( (((n_pts+1)/2) - ((bin_count-1)/2)) : (((n_pts+1)/2) + ((bin_count-1)/2)));
                end
            end
            x_vals(bin==k) = spacings;
        end
        all_x{i}(:,j) = x_vals/2.5*swarm_width + base_spacing(idx);
    end
end

% PLOT datapoints

% Create figure (if axes wasn't provided)
if ~ishandle(p.Results.Axes)
    swarmfig = figure('Position', positionfig(600,400), 'PaperPositionMode','auto');
    swarmaxes = axes('Parent',swarmfig);
else
    swarmaxes = p.Results.Axes;
end

% Place points/ bars in figure
hold(swarmaxes,'on')
for i = 1:length(data)
    for j = 1:size(data{i},2)
        if n_subgroups==1
            clr1 = clr(i,:);
        else
            clr1 = clr(j,:);
        end
        
        plot(all_x{i}(:,j),data{i}(:,j),'o','MarkerFaceColor',clr1,'MarkerEdgeColor',clr1/2,...
            'MarkerSize',marker_size)
    end
end


% Place means (if set)
if show_means
    idx = 0;
    for i = 1:length(data)
        for j = 1:size(data{i},2)
            idx = idx+1;
            plot([base_spacing(idx)-0.3*swarm_width, base_spacing(idx)+0.3*swarm_width], [mean(data{i}(:,j)), mean(data{i}(:,j))],...
                    'Color',[0.4431 0.4510 0.4627],'LineWidth',3)
        end
    end
end

hold(swarmaxes,'off')




set(swarmaxes,'XLim',[min(base_spacing)-1,max(base_spacing)+1],'XTick',[])




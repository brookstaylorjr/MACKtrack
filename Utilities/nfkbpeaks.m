function [peak_times, peak_amps, valley_times, valley_amps, fig1] = nfkbpeaks(trajectories, varargin)
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [peak_times, peak_amps, valley_times, valley_amps, fig1] = nfkbpeaks(trajectories, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBPEAKS uses globalpeaks, plus some additional filtering steps, to find peaks in a trajectory of NFkB activity
%
% INPUTS (required)
% trajectories     matrix of cell trajectories (each row is a cell)
%
% INPUTS (optional - name/value pairs)
% 'ShowGraph'      boolean flag to show sampling of found peaks on graphed trajectories (2x6 examples) - default is 0
% 'NumPeaks'       peaks to look for in trajectory - default is 16
% 'BeginFrame'     supress peaks that are found before this point (default =3)
% 'EndFrame'       supress peaks that are found after this point (default = cols in trajectories - 2)
% 'MinDist'        supress peaks that are within n frames of one another (default = 6)
% 'MinHeight'      supress peaks whose absolute amplitude falls below this value (default = 0.75)
% 'SmoothSize'     window to smooth trajectories slightly (using 'smoothrows' - default = 3 frames)
% 'YLim'           y limits for graph (default = [0.25 9])

% OUTPUTS
% peak_times       locations of all peaks (column indicies)
% peak_amps        heights of all peaks
% valley_times     locations of all valleys (column indicies)
% valley_amps      heights of all valleys
% fig1             handle of (optional) output figure
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% Create input parser object, add required params from function input
p = inputParser;
% Required: graph_data
valid_data = @(x) assert((isnumeric(x)&&ismatrix(x)),...
    'input ''trajectories'' must be a matrix (where rows correspond to a set measurements for one individual');
addRequired(p,'trajectories',valid_data);

% Optional parameters
valid_val = @(x) assert(isnumeric(x)&&(numel(x)==1), 'input must be single numeric value');
valid_lim = @(x) assert(isnumeric(x)&&(numel(x)==2), 'input must be vector specifying [y_min y_max]');

addParameter(p,'ShowGraph',0,valid_val);
addParameter(p,'NumPeaks',16,valid_val);
addParameter(p,'BeginFrame',3,valid_val);
addParameter(p,'EndFrame',size(trajectories,2)-2,valid_val);
addParameter(p,'MinHeight',0.75,valid_val);
addParameter(p,'MinDist',9,valid_val);
addParameter(p,'SmoothSize',3,valid_val);
addParameter(p,'YLim',[0.25 9] ,valid_lim);



% Parse parameters, assign to variables
parse(p,trajectories, varargin{:})
num_peaks = p.Results.NumPeaks;
begin_frame = p.Results.BeginFrame;
end_frame = p.Results.EndFrame;
min_dist = p.Results.MinDist;
min_height = p.Results.MinHeight;


%% Parameters (peak attributes)

% Initialize empty data (per cell)
peak_times = nan(size(trajectories,1),num_peaks);
peak_amps = nan(size(trajectories,1),num_peaks);
valley_times = nan(size(trajectories,1),num_peaks-1);
valley_amps = nan(size(trajectories,1),num_peaks-1);


% Sample graph (show peak-finding performace)
if p.Results.ShowGraph
    fig1 = figure('PaperPositionMode','auto','Position',positionfig(920, 94));
    ha = tight_subplot(2,6,[0.01 0.01]);
end

% Smooth data slightly
if p.Results.SmoothSize>1
    nfkb_smooth = smoothrows(trajectories,p.Results.SmoothSize);
else
    nfkb_smooth = trajectories;
end

% Cycle through each cell; find peaks and plot (if selected)
plot_idx=1;
for i = 1:size(peak_times,1)
    vect = nfkb_smooth(i,1:min([end,end_frame]));
    if sum(~isnan(vect)) > 10
        [pks, locs, heights] = globalpeaks(vect,num_peaks);
    else
        pks = []; locs = []; heights = [];
    end
    % Supress any peaks that are too close to one another
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    heights = heights(order);
    while min(diff(locs)) < min_dist
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (heights(tmp)>=heights(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];
        heights(tmp) = [];
    end
    
    % Supress early and low-quality peaks
    drops = (locs<begin_frame) | (heights<min_height);
    pks(drops) = [];
    locs(drops) = [];
    if ~isempty(locs)
        peak_times(i,1:length(locs)) = locs;
        peak_amps(i,1:length(locs)) = pks;
        % Add 1st two valleys between consecutive peaks
        for j = 1:(length(locs)-1)
            valley_loc = find(vect(locs(j):locs(j+1))==min(vect(locs(j):locs(j+1))),1,'first');
            valley_times(i,j) = valley_loc+locs(j)-1;
            valley_amps(i,j) = vect(valley_loc+locs(j)-1);
        end
    end
    % Show small multiples with called peaks - some randomness is introduced so same cells aren't always shown.
    % (if graph isn't filling up, increase cutoff (0.2 seems to work ok for 300-500 cells)
    if p.Results.ShowGraph
        if (plot_idx<=length(ha)) && (length(pks)>6) && (rand(1) < 0.2)
            plot(ha(plot_idx),1:length(vect),vect)
            hold(ha(plot_idx),'on')
            plot(ha(plot_idx),locs,pks,'ok','MarkerSize',4)
            %plot(ha(plot_idx),valley_times(i,:),valley_amps(i,:),'or','MarkerSize',4)
            hold(ha(plot_idx),'off')
            set(ha(plot_idx),'XTickLabel',{},'YTickLabel',{},'XLim',[1 end_frame],'YLim',p.Results.YLim)
            plot_idx = plot_idx+1;
        end
    end
end

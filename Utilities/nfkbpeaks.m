function [peak_times, peak_amps, valley_times, valley_amps, fig1] = nfkbpeaks(cell_trajectories, show_graph)
%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [peak_data, valley_data] = nfkbpeaks(cell_trajectories, show_graph)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBPEAKS uses globalpeaks, plus some additional filtering steps, to find peaks in a trajectory of NFkB activity
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default: don't no graph output
if nargin<2
    show_graph = 0;
end

%% Parameters (peak attributes)
num_peaks = 16; % Maximum number of peaks we can find from each trajectory
begin_frame = 3; % Start looking for peaks @ this frame
end_frame = 216; % Stop looking for peaks @ this frame
min_dist = 9; % Minimum allowable distance between peaks
min_height = 0.75; % Minimum height to be considered a peak
smooth_size = 3; % Size of sliding window to smooth trajectories slightly
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Initialize empty data (per cell)
peak_times = nan(size(cell_trajectories,1),num_peaks);
peak_amps = nan(size(cell_trajectories,1),num_peaks);
valley_times = nan(size(cell_trajectories,1),2);
valley_amps = nan(size(cell_trajectories,1),2);

% Sample graph (show peak-finding performace)
if show_graph
    fig1 = figure('PaperPositionMode','auto','Position',positionfig(880, 90));
    ha = tight_subplot(2,6,[0.01 0.01]);
end

% Smooth data slightly
nfkb_smooth = smoothrows(cell_trajectories,smooth_size);

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
        if length(locs)>2
            valley_loc = find(vect(locs(1):locs(2))==min(vect(locs(1):locs(2))),1,'first');
            valley_times(i,1) = valley_loc+locs(1)-1;
            valley_amps(i,1) = vect(valley_loc+locs(1)-1);
            if length(locs)>3
                valley_loc = find(vect(locs(2):locs(3))==min(vect(locs(2):locs(3))),1,'first');
                valley_times(i,2) = valley_loc+locs(2)-1;
                valley_amps(i,2) = vect(valley_loc+locs(2)-1);
            end
        end
    end
    % Show small multiples with called peaks - some randomness is introduced so same cells aren't always shown.
    % (if graph isn't filling up, increase cutoff (0.2 seems to work ok for 300-500 cells)
    if show_graph
        if (plot_idx<=length(ha)) && (length(pks)>6) && (rand(1) < 0.2)
            plot(ha(plot_idx),1:length(vect),vect)
            hold(ha(plot_idx),'on')
            plot(ha(plot_idx),locs,pks,'ok','MarkerSize',4)
            hold(ha(plot_idx),'off')
            set(ha(plot_idx),'XTickLabel',{},'YTickLabel',{},'XLim',[1 end_frame],'YLim',[-0.25 9])
            plot_idx = plot_idx+1;
        end
    end
end

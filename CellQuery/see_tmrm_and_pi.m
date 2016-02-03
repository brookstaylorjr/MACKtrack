function [graph, info, measure] = see_tmrm_and_pi(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_tmrm_and_pi(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_TMRM is a visualization function to see measured TMRM dye 
%
%
% id             experiment ID (from Google Spreadsheet specigied in "loadID.m")
% show_graphs    boolean flag; specifies whether standard behavioral graphs will be shown
% diagnos        boolean flag; specifies whether optional diagnostic graphs will be shown
%
% 
% graph          primary output structure; must specify
%                   1) filtered/processed data (graph.var) 
%                   2) time vector for all images (graph.t) 
% info           secondary output structure; must specify
%                   1) Y limits for graphing (info.graph_limits)
%                   2) parameters from loadID.m (info.parameters) 
% measure         full output structure from loadID
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


%% Setup
if nargin<3
    diagnos=0;
    if nargin<2
        show_graphs = 0;
    end
end
%%

% Load data; set parameters
[measure, info] =loadID(id);
info.parameters.FramesPerHour = 40; % 1.5 min between frames
info.Module = 'tmrmModule';
t_max = (size(measure.TMRM_cytoMean,2)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
info.graph_limits = [-20 120];
info.graph_limits2 = [200 400];



%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(measure.TMRM_cytoMean(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.TMRM_cytoMean(:,1:120)),2)>3]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;


%% Outputs
% Extract TMRM measurement, apply filtering, and subtract baseline for each cell
graph.var = measure.TMRM_cytoMean(info.keep,:);
baseline = nanmean(graph.var(:,1:4),2); % caluclate baseline from avg of 1st 4 frames
graph.var = graph.var - repmat(baseline,1, size(graph.var,2)); % copy baseline for each timepoint and subtract

% Extract PI data (no normalization)
graph.var2 = measure.MeanIntensityNuc(info.keep,:);


graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    figure,imagesc(graph.var,info.graph_limits)
    figure,imagesc(graph.var2,info.graph_limits2)
end
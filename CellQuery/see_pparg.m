function [graph, info, measure] = see_pparg(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_pparg(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PPARG is a visualization function to see measured nuclear PPARg levels in tagged cells
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


% Load data; set parameters
[measure, info] = loadID(id);
info.parameters.FramesPerHour = 6; % 10 min between frames
info.Module = 'nucintensityModule';
t_max = (size(measure.MeanIntensityNuc,2)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
info.graph_limits = [400 1000];



%% Filtering
droprows = zeros(size(measure.MeanIntensityNuc,1),1);
%droprows = [droprows, sum(isnan(measure.MeanIntensityNuc(:,1:4)),2)>2]; % Cells existing @ expt start
%droprows = [droprows, sum(isnan(measure.MeanIntensityNuc(:,1:120)),2)>3]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;


%% Outputs
% Extract measurement and apply filtering

all_pparg = measure.MeanIntensityNuc(info.keep,:);
% Normalize by image background (per position)
for i = 1:length(info.parameters.XYRange)
    xy = info.parameters.XYRange(i);
    all_pparg(info.CellData(:,1)==xy,:) = (all_pparg(info.CellData(:,1)==xy,:) - info.parameters.adj_distr(1,i)) ;
    
end


graph.var = all_pparg;

graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    figure,imagesc(graph.var,info.graph_limits)
end
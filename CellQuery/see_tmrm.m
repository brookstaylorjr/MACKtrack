function [graph, info, measure] = see_tmrm(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_tmrm(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_TMRM is a visualization function to measure the TMRM dye. 
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
%                   3) XY convection adjustment (graph.shift) 
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
info.parameters.FramesPerHour = 40; % 1.5 min between frames
info.Module = 'intensityModule';
t_min = (size(measure.TMRM_cytoMean,2)-1)/info.parameters.FramesPerHour*60; % Number of hours to display in graphs
info.graph_limits = [400 750];

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(measure.Area(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.TMRM_cytoMean(:,1:60)),2)>3]; % Long-lived cells



%% Outputs
info.keep = max(droprows,[],2) == 0;
graph.var = measure.TMRM_cytoMean(info.keep,:);
graph.t = 0:(60/info.parameters.FramesPerHour):t_min;
graph.var = graph.var(:,1:length(graph.t));


graph.celldata = info.CellData(info.keep,:);
graph.opt = maketicks(graph.t,info.graph_limits,0);
graph.opt.Name = 'Number of identified spots'; 
graph.shift = zeros(length(unique(info.CellData(:,1))),1);
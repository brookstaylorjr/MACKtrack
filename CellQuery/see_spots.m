function [graph, info, measure] = see_spots(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_spots(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_SPOTS is a visualization function for 'NumSpots' in 'SpotCountModule'. 
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
info.ImageExpr = info.parameters.spotCountModule.ImageExpr;
t_hrs = min([21,size(measure.NumSpots,2)/info.parameters.FramesPerHour]); % Number of hours to display in graphs
info.graph_limits = [0 20];

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(measure.Area(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.NumSpots(:,1:100)),2)>3]; % Long-lived cells



%% Outputs
info.keep = max(droprows,[],2) == 0;
graph.var = measure.NumSpots(info.keep,:);
graph.t = 0:(1/info.parameters.FramesPerHour):t_hrs;
graph.var = graph.var(:,1:length(graph.t));


graph.celldata = info.CellData(info.keep,:);
graph.opt = maketicks(graph.t,info.graph_limits,0);
graph.opt.Name = 'Number of identified spots'; 
graph.shift = zeros(length(unique(info.CellData(:,1))),1);
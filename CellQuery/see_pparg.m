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
info.Module = 'ppargModule';
t_max = (size(measure.MedianPPARg,2)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs


all_pparg = measure.MeanPPARg;
info.graph_limits = prctile(all_pparg(~isnan(all_pparg)),[3 97]);



% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_pparg;
find_parent = @(row) find((info.CellData(:,1) == row(1)) & (info.CellData(:,2)== row(5)));
graph.divisions = [];
graph.divide_pts = false(size(graph.var_all));
for i = 1:size(graph.var_all,1)
    if info.CellData(i,5)>0      
        graph.var_all(i,1:info.CellData(i,3)) = graph.var_all(find_parent(info.CellData(i,:)),1:info.CellData(i,3));
        graph.divide_pts(find_parent(info.CellData(i,:)),info.CellData(i,3)) = 1;
        graph.divide_pts(i,1:info.CellData(i,3)) = graph.divide_pts(find_parent(info.CellData(i,:)),1:info.CellData(i,3));
        graph.divisions = cat(1,graph.divisions,[i, info.CellData(i,3)]);
    end
end
all_pparg = graph.var_all;

%% Filtering
droprows = zeros(size(all_pparg,1),1);
droprows = [droprows, sum(isnan(all_pparg(:,end-3:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_pparg),2)>50]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;


%% Outputs
% Extract measurement and apply filtering
all_pparg = all_pparg(info.keep,:);
graph.divide_pts = graph.divide_pts(info.keep,:);

graph.var = all_pparg;

graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    figure,imagesc(graph.var,info.graph_limits)
end
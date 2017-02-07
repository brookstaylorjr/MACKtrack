function [graph, info, measure] = see_Hoechst(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_pi(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PI is a visualization function to see measured vital cell dye levels
% (like P.I. or sytox green)
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
info.parameters.FramesPerHour = 40; % 1.5 min between frames
info.Module = 'nucintensityModule';
t_max = (size(measure.MeanIntensityNuc2,2)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
info.graph_limits = [200 400];



%% Filtering
droprows = [];
%include everything!
%droprows=zeros(size(measure.MeanIntensityNuc2,1),1);
droprows = [info.droprows, info.CellData(:,end)]; % Filter out edge cells
%droprows = [droprows, sum(isnan(measure.MeanIntensityNuc2(:,1:4)),2)>2]; % Cells existing @ expt start
%droprows = [droprows, sum(isnan(measure.MeanIntensityNuc2(:,1:120)),2)>3]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;

%% Outputs
% Extract measurement and apply filtering
graph.var = measure.MeanIntensityNuc2(info.keep,:);
graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    figure,imagesc(graph.var,info.graph_limits)
end
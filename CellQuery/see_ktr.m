function [graph, info, measure] = see_ktr(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_ktr(id,show_graphs, diagnos)
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
info.Module = 'ktrModule';
graph.t = (measure.MultiKTR_t(1,:)-1)/info.parameters.FramesPerHour; % Number of hours to display in graphs


all_ktr = measure.KTR_ratio1;
all_ktr(all_ktr==0) = nan;


%all_ktr = 1./all_ktr;

% Add parent trajectories to children
all_ktr = copychildren(all_ktr,info.CellData, measure.MultiKTR_t(1,:));



%% Filtering
droprows = zeros(size(all_ktr,1),1);
droprows = [droprows, sum(isnan(all_ktr(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(all_ktr(:,1:300)),2)>100]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;


%% Outputs
% Extract measurement and apply filtering
all_ktr = all_ktr(info.keep,:);

%all_ktr = all_ktr - repmat(nanmean(all_ktr(:,1:4),2),[1,size(all_ktr,2)]);
info.graph_limits = prctile(all_ktr(~isnan(all_ktr)),[10 90]);

graph.var = all_ktr;


graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    colormaps=loadcolormaps;
    figure,imagesc( [min(graph.t) max(graph.t)],[1 size(graph.var,1)], graph.var)
    colormap(colormaps.magma)
    set(gca,'Clim',info.graph_limits)
end
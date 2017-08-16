function [graph, info, measure] = see_fret(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_pparg(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_FRET is a visualization function to see measured cytoplasmic FRET levels in tagged cells
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
[measure, info] = loadID(id);
info.ImageExpr = ['makefret(',info.parameters.fretModule.ImageExpr,')'];

t_max = (size(measure.MeanFRET_cyto,2)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
info.graph_limits = [0.7 2];
all_fret = measure.MeanFRET_cyto;

% Add parent trajectories to children
find_parent = @(row) find((info.CellData(:,1) == row(1)) & (info.CellData(:,2)== row(5)));
for i = 1:size(all_fret,1)
    if info.CellData(i,5)>0
        all_fret(i,1:info.CellData(i,3)) = all_fret(find_parent(info.CellData(i,:)),1:info.CellData(i,3));
    end
end

graph.var_all = all_fret;
[all_fret, graph.lineage] = copychildren(all_fret, info.CellData);



%% Filtering
droprows = zeros(size(measure.MeanFRET_cyto,1),1);
%droprows = [droprows, sum(isnan(measure.MeanIntensityNuc(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(all_fret(:,1:min([200 size(all_fret,2)]))),2)>15]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;

%% Subtract basal level of FRET (?)

all_fret = all_fret(info.keep,:);

all_fret_baseline  = all_fret - repmat(nanmean(all_fret(:,1:2),2),[1 size(all_fret,2)]);



%% Outputs
% Extract measurement and apply filtering
graph.var = all_fret;
graph.var2 = all_fret_baseline;

graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

graph.celldata = info.CellData(info.keep,:);
graph.shift = zeros(length(unique(info.CellData(:,1))),1);

if show_graphs
    figure,imagesc(graph.var,info.graph_limits)
end
function [graph, info, measure] = see_paprg(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_paprg(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PPARG is a basic visualization function to plot single-cell expression levels of PPARg (or similar) over time
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'         'on' or 'off' - show verbose output
%
% OUTPUTS:  
% graph          primary output structure; must specify
%                   1) filtered/processed data (graph.var) 
%                   2) time vector for all images (graph.t) 
%                   3) XY convection adjustment (graph.shift) 
% info           secondary output structure; must specify
%                   1) Y limits for graphing (info.graph_limits)
%                   2) parameters from loadID.m (info.parameters) 
% measure         full output structure from loadID
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||isstruct(x)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'MinLifetime',100, @isnumeric);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end

% Load data; set parameters
[measure, info] = loadID(id);
t_max = (length(info.parameters.TimeRange)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
if isfield(measure,'MeanPPARg')
    all_pparg = measure.MeanPPARg;
    info.ImageExpr = info.parameters.ppargModule.ImageExpr;
else
    all_pparg = measure.MeanNuc1;
    info.ImageExpr = info.parameters.intensityModule.ImageExpr;
end

info.graph_limits = prctile(all_pparg(~isnan(all_pparg)),[3 97]);


% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_pparg;
[all_pparg, graph.lineage] = copychildren(all_pparg, info.CellData);

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(all_pparg(:,end-3:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_pparg),2)>50]; % Long-lived cells
info.keep = max(droprows,[],2) == 0;
% Show some filter information

if verbose_flag
    filter_str = {'didn''t exist @ experiment end', 'short-lived cells'};
    disp(['INITIAL: ', num2str(size(droprows,1)),' cells'])
    
    for i = 1:size(droprows,2)
        if i ==1
            num_dropped = sum(droprows(:,i)==1);
        else
            num_dropped = sum( (max(droprows(:,1:i-1),[],2)==0) & (droprows(:,i)==1));
        end
        disp(['Filter #', num2str(i), ' (',filter_str{i},') - ',num2str(num_dropped), ' cells dropped']) 
    end
    disp(['FINAL: ', num2str(sum(max(droprows,[],2) == 0)),' cells'])
end

%% Outputs -> drop filtered cells 
graph.var = all_pparg(info.keep,:);
graph.lineage = graph.lineage(info.keep,:);
graph.celldata = info.CellData(info.keep,:);
graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

if graph_flag
    figure,imagesc(graph.var,info.graph_limits)
end
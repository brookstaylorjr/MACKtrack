function [graph, info, measure] = see_h2b(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_h2b(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_H2B is a basic visualization function to plot single-cell expression levels of H2B (or similar) over time.
% If a 'Measurement' is not provided, this function will default to 'IntegratedNuc3'.
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Measurement'     name of measurement field that will be read as PPARg - defaults to 'MeanNuc1'
% 'ImageExpr'       expr. for the fluorescence image used to measure PPARg - defaults to 1st slot of intensityModule
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'         'on' or 'off' - show verbose output
%
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
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||isstruct(x)||exist(x,'file')||exist(x,'dir'),...
    'ID input must be either: (1) an AllMeasurements structure, (2) AllMeasurements file or folder location, or (3) a spreadsheet ID');
addRequired(p,'id',valid_id);

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Measurement','', @ischar);
addParameter(p,'ImageExpr','', @ischar);
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end
measure_field = p.Results.Measurement;

% Load data; set parameters
[measure, info] = loadID(id);
t_max = (length(info.parameters.TimeRange)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs

if ~isempty(measure_field)
    assert(isfield(measure,measure_field),['Error: supplied measurement field (',measure_field,') not found.'])
    assert(~isempty(p.Results.ImageExpr),['A matching image expression must be provided for measurement ', measure_field])
    info.ImageExpr = p.Results.ImageExpr;
else
    if isfield(measure,'IntegratedNuc3')
        measure_field = 'IntegratedNuc3';
        info.ImageExpr = info.parameters.intensityModule.ImageExpr3;
    else
        error('No valid H2B measurement found')
    end
end
all_h2b = measure.(measure_field);
info.graph_limits = prctile(all_h2b(~isnan(all_h2b)),[3 97]);


% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_h2b;
[all_h2b, graph.lineage] = copychildren(all_h2b, info.CellData);

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(all_h2b(:,end-3:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_h2b),2)>50]; % Long-lived cells
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
graph.var = all_h2b(info.keep,:);
graph.lineage = graph.lineage(info.keep,:);
graph.celldata = info.CellData(info.keep,:);
graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

if graph_flag
    figure,imagesc(graph.var,info.graph_limits)
end
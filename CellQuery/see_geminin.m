function [graph, info, measure] = see_geminin(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_geminin(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PPARG is a basic visualization function to plot single-cell
% expression levels of geminin
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Measurement'     name of measurement field that will be read as PPARg - defaults to 'MeanNuc1'
% 'ImageExpr'       expr. for the fluorescence image used to measure PPARg - defaults to 1st slot of intensityModule
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'         'on' or 'off' - show verbose output
% 'MinLifetime'     cell trajectories with < MinLifetime frames will be filtered out. Default = 100.
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
addParameter(p,'Measurement','MeanNuc2', @ischar);
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'MinLifetime',200, @isnumeric);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end
measure_field = p.Results.Measurement;

% Load data; set parameters
[measure, info] = loadID(id);
t_max = (length(info.parameters.TimeRange)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs

all_geminin = measure.(measure_field);
all_pparg = measure.MeanNuc1;

info.ImageExpr = info.parameters.intensityModule.ImageExpr;




% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_geminin;
[all_geminin, graph.lineage] = copychildren(all_geminin, info.CellData);

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(all_geminin(:,end-5:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_geminin),2)>p.Results.MinLifetime]; % Long-lived cells
droprows = [droprows, sum(isnan(all_geminin),2) > 2];

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


%% Process data
all_geminin = all_geminin - repmat(prctile(all_geminin(:,1:100),4,2),[1 size(all_geminin,2)]);

%% Outputs -> drop filtered cells

info.graph_limits = prctile(all_geminin(~isnan(all_geminin)),[3 97]);


graph.var = all_geminin(info.keep,:);
graph.lineage = graph.lineage(info.keep,:);
graph.celldata = info.CellData(info.keep,:);
graph.t = 0:(60/info.parameters.FramesPerHour):t_max;

if graph_flag
    figure,imagesc(graph.var,info.graph_limits)
end
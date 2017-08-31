function [graph, info, filtered_cells, measure] = see_paprg_EBM(id,Mt_window,alpha_percentile,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_paprg_EBM(id,Mt_window,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PPARG is a basic visualization function to plot single-cell expression levels of PPARg (or similar) over time
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
% Mt_window      number of frames analyzed together by MACKtrack
% alpha_percentile        percentile to cutoff probable mistraces
%
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
% filtered_cells numbers of cells filtered out at each filtering step                 
% measure         full output structure from loadID
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||isstruct(x)||exist(x,'file')||exist(x,'dir'),...
    'ID input must be either: (1) an AllMeasurements structure, (2) AllMeasurements file or folder location, or (3) a spreadsheet ID');
addRequired(p,'id',valid_id);

% Required: Mt_window
if exist('Mt_window') == 0
    Mt_window = 13;
end

% Required: alpha_percentile
if exist('alpha_percentile') == 0
    alpha_percentile = 95;
end

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Measurement','MeanNuc1', @ischar);
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'MinLifetime',100, @isnumeric);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end
measure_field = p.Results.Measurement;

% Load data; set parameters
[measure, info] = loadID(id);
t_max = (length(info.parameters.TimeRange)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs
if isfield(measure,'MeanPPARg')
    all_pparg = measure.MeanPPARg;
    info.ImageExpr = info.parameters.ppargModule.ImageExpr;
else
    all_pparg = measure.(measure_field);
    info.ImageExpr = info.parameters.intensityModule.ImageExpr;
end

info.graph_limits = prctile(all_pparg(~isnan(all_pparg)),[3 97]);

% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_pparg;
[all_pparg, graph.lineage] = copychildren(all_pparg, info.CellData);

%% Compute variables for the filtering of irreversible PPARg jumps
% Integral of the difference between PPARg and the local PPARg average
% squared
[rows, cols] = size(all_pparg);
integral = zeros([rows, cols]);
for k = Mt_window:(cols-Mt_window)
    for l = 1:rows
        local_integral = 0;
        local_average = 0;
        for m = 1:(2*Mt_window)
            local_average = local_average + all_pparg(l,k+m-Mt_window)/(2*Mt_window);
        end
        for n = 1:(2*Mt_window)
            local_integral = local_integral + (all_pparg(l,k+n-Mt_window)-local_average)^2;
        end
        integral(l,k) = local_integral/(2*Mt_window);
    end
end
max_integral = max(integral,[],2);
alpha = prctile(max_integral,alpha_percentile);

%Show a histogram of all integral values to chose cut-off (alpha)
% figure
% histogram(max_integral);
% hold on
% xlabel('Alpha');
% ylabel('# cells');
% xlim([0 1000]);

%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(all_pparg(:,end-3:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_pparg),2)>10]; % Long-lived cells
droprows = [droprows, (max_integral>alpha)>0]; %PPARg jumps
info.keep = max(droprows,[],2) == 0;
% Show some filter information

if verbose_flag
    filter_str = {'didn''t exist @ experiment end', 'short-lived cells','mistraces based on PPARg jumps'};
    disp(['INITIAL: ', num2str(size(droprows,1)),' cells'])

filtered_cells = zeros([1 size(droprows,2)]);
    for i = 1:size(droprows,2)
        if i ==1
            num_dropped = sum(droprows(:,i)==1);
        else
            num_dropped = sum( (max(droprows(:,1:i-1),[],2)==0) & (droprows(:,i)==1));
        end
        filtered_cells(i) = num_dropped;
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
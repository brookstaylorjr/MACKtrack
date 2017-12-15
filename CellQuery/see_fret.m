function [graph, info, measure] = see_fret(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_fret(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_PPARG is a basic visualization function to plot single-cell expression levels of PPARg (or similar) over time.
% If a 'Measurement' is not provided, this function will use 'MeanPPARg' or 'MeanNuc1'.
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'         'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'         'on' or 'off' - show verbose output
% 'GetImage'        if provided here, see_fret will instead return a FRET image corresponding to [site timept]
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
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'GetImage',[nan nan], @isnumeric);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end

% Load data; set parameters
[measure, info] = loadID(id,0);
t_max = (length(info.parameters.TimeRange)-1)/(info.parameters.FramesPerHour/60); % Number of hours to display in graphs


%% Alternate function mode: use either FRET module or multi_FRET module to return an appropriate image
if sum(isnan(p.Results.GetImage)) == 0
    if info.parameters.fretModule.Use
        expr1 = info.parameters.fretModule.ImageExpr;
        expr2 = info.parameters.fretModule.ImageExpr2;
    elseif info.parameters.multi_fretModule.Use
        expr1 = info.parameters.multi_fretModule.ImageExpr;
        expr2 = info.parameters.multi_fretModule.ImageExpr2;
    elseif info.parameters.fret_lowModule.Use
        expr1 = info.parameters.fret_lowModule.ImageExpr;
        expr2 = info.parameters.fret_lowModule.ImageExpr2;
    end
    
    
    % Image correction: use rough version of flatfield correction (divide by raw FF image
    i = p.Results.GetImage(1);    j = round(100*p.Results.GetImage(2))/100;
    fret = checkread(namecheck([info.locations.scope,filesep,info.parameters.ImagePath,filesep,eval(expr1)]));
    fret = flatfieldcorrect(fret,double(info.parameters.Flatfield{1}));    
    fret = fret- prctile(fret(:),2);
    
    cfp = checkread(namecheck([info.locations.scope,filesep,info.parameters.ImagePath,filesep,eval(expr2)]));
    cfp = flatfieldcorrect(cfp,double(info.parameters.Flatfield{1}));  
    cfp = cfp - prctile(cfp(:),2);
    graph = fret./cfp;
    
    graph(cfp<16) = 0.1;
    graph(graph<0 | graph>10) = 0.1;
    
    info = [];
    measure = [];
    return;
end


%%
if ~strcmpi(info.parameters.ImageType,'None')
    all_fret = measure.FRETpctile_cyto(:,:,7);
else
    all_fret = measure.MedianFRET_nuc;
end

info.graph_limits = prctile(all_fret(~isnan(all_fret)),[3 97]);

info.ImageExpr = '';


% Add parent trajectories to children - record time/index of divisions.
graph.var_all = all_fret;
[all_fret, graph.lineage] = copychildren(all_fret, info.CellData);


%% Filtering
droprows = [];
droprows = [droprows, sum(isnan(all_fret(:,1:5)),2)>2]; % Cells existing @ expt beginning
droprows = [droprows, sum(isnan(all_fret(:,end-5:end)),2)>2]; % Cells existing @ expt end
droprows = [droprows, sum(isnan(all_fret),2)>50]; % Make sure there are not many missing vals
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
graph.var = all_fret(info.keep,:);
graph.lineage = graph.lineage(info.keep,:);
graph.celldata = info.CellData(info.keep,:);

if ~isfield(measure,'MultiFRET_t')
    graph.t = 0:(60/info.parameters.FramesPerHour):t_max;
else
    graph.t = (measure.MultiFRET_t-1)/info.parameters.FramesPerHour;
    % If I want to see all timepoints, need to back-convert timept into (e.g.) 100.01, etc.
    t = measure.MultiFRET_t;
    z_idx = [find(mod(t,1)==0),length(t)+1];
    graph.t_aux = zeros(size(t));
    for i = 1:(length(z_idx)-1)
        info.t_aux(z_idx(i):z_idx(i+1)-1) = t(z_idx(i)) + ((1:(z_idx(i+1)-z_idx(i)))-1)/100;
    end
end

if graph_flag
    figure,imagesc(graph.var,info.graph_limits)
end


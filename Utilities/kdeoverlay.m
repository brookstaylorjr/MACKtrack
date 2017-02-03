function [h_ax, bandwidth] = kdeoverlay(vects, varargin) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = kdeoverlay(vects, places, varargin) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% KDEOVERLAY computes a kernel density estimate for a set of distributions ('vects')
%
% INPUTS (required)
% vects          1xN cell array of vectors (1-D array of object measurements)
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Color'        1x... cell vector specifying area fill colors - cycles if length < N
% 'Xlim'         2 element vector with graph [x_min, x_max]. Default is 1st and 99th percentile of all data
% 'Bandwidth'    Bandwidth of kernel density estimate (if not provided, will use MATLAB default)
% 'Axes'         Axes handle of axes where new violin figure will be plotted (default: create new figure)
% 'LineWidth'    Line width around shape - default = 1 (can be 0)
% 'Alpha'        Transparency of shape(s). 1 = fully opaque. 0 = fully transparent

%
% OUTPUTS
% violin        Axes handle of violin figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% 1) Vector data (must be cell matric)
addRequired(p,'vects',@iscell);

% Optional parameters
default_color = mat2cell(linspecer(numel(vects)),ones(1,numel(vects)));
valid_color = @(x) assert(iscell(x)&&length(x{1})==3, 'Specify colors with a cell matrix of RGB triplets');
addParameter(p,'Color', default_color,valid_color);

all = cell2mat(vects(:));
valid_xlim = @(x) assert(length(x)==2,'XLim must be a 2 element vector');
addParameter(p,'XLim', prctile(all(:),[1 99]),valid_xlim);
addParameter(p,'Bandwidth',nan,@isnumeric);
addParameter(p,'Axes',nan,@ishandle);
valid_width = @(x) assert(isnumeric(x)&&(x>=0),'Line width must be >= 0');
addParameter(p,'LineWidth',1,valid_width);
valid_alpha = @(x) assert(isnumeric(x)&&(x>=0)&&(x<=1),'Alpha must be between 0.0 and 1.0 (inclusive)');
addParameter(p,'Alpha',0.4,valid_width);


% Parse inputs, save some to variables
parse(p,vects, varargin{:})
xlim = p.Results.XLim;
colors = p.Results.Color;
bandwidth = p.Results.Bandwidth;
linewidth = p.Results.LineWidth;
alpha = p.Results.Alpha;

% Create figure (if axes wasn't provided)
if ~ishandle(p.Results.Axes)
    kdefig = figure('Position', [500, 1031, 800, 300], 'PaperPositionMode','auto');
    h_ax = axes('Parent',kdefig);
else
    h_ax = p.Results.Axes;
end
%%
% If bandwidth isn't defined, get an estimate
x_val = linspace(min([xlim(:);all]) ,max([xlim(:);all]),1024);
if isnan(bandwidth)
    [~,~,bandwidth] = ksdensity(all,x_val,'function','pdf');
end

hold(h_ax,'on')
for i = 1:length(vects)
    [y_val] = ksdensity(vects{i},x_val,'function','pdf','bandwidth',bandwidth);
    
    if linewidth>0
        fill([x_val,fliplr(x_val)],[y_val,zeros(size(y_val))],colors{i},'FaceAlpha',alpha,'EdgeColor',colors{i},'LineWidth',linewidth)
    else
        fill([x_val,fliplr(x_val)],[y_val,zeros(size(y_val))],colors{i},'FaceAlpha',alpha,'EdgeColor','none')
    end
end

set(h_ax,'XLim',xlim)






function violin = violin(vects, places, varargin) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = violin(vects, places, varargin) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% VIOLIN creates a violin plot, spacing them according to a secondary vector (e.g. doses)
%
% INPUTS (required)
% vects          1xN cell array of vectors (1-D array of object measurements)
% places         1xN array directing placement of each Violin
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Color'        1x... cell vector specifying violin fill colors - cycles if length < N
% 'YLim'         2 element vector with graph [y_min, y_max]. Default is 5th and 95th percentile of all data
% 'ShowBins'     Show additional histogram  graph ('on' or 'off' - default is 'off')
% 'Area'         Total area of graph taken up by each shape (default = 0.01)
% 'BinScale'     Scaling factor (from default) for number of histogram bins (can specify scalar or a vector of length N)
% 'Bins'         Vector of bin centers - if provided, 'BinScale' will be ignored
% 'XSpace'       Axis whitespace before the first violin plot, and after the last (default = 0.1)
% 'Axes'         Axes handle of axes where new violin figure will be plotted (default: create new figure)
% 'Smoothing'    Smoothing of violin shapes (such that histogram bars are evident, or smoothed out). Default = 'on'
% 'Connect'      Add connecting line between shapes (default = 'on')
% 'LineWidth'    Line width around shape - default = 0.5 (can be 0)
%
% OUTPUTS
% violin        Axes handle of violin figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% 1) Vector data (must be cell matric)
addRequired(p,'vects',@iscell);
% 2) X axis Placement
valid_places = @(x) assert(numel(x)==numel(vects), '2nd argument must be same size as 1st');
addRequired(p,'places',valid_places);

% Optional parameters
colors = setcolors;
default_color = colors.peacock(end:-1:1);

valid_color = @(x) assert(iscell(x)&&length(x{1})==3, 'Specify colors with a cell matrix of RGB triplets');
addParameter(p,'Color', default_color,valid_color);

all = cell2mat(vects(:));
valid_ylim = @(x) assert(length(x)==2,'YLim must be a 2 element vector');
addParameter(p,'YLim', prctile(all(:),[5 95]),valid_ylim);
expectedFlags = {'on','off'};
addParameter(p,'ShowBins','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Area',0.01,@isnumeric);
addParameter(p,'XSpace',0.1,@isnumeric);
addParameter(p,'BinScale',1,@isnumeric);
valid_width = @(x) assert(isnumeric(x)&&(x>=0),'Line width must be >= 0');
addParameter(p,'LineWidth',0.5,valid_width);
valid_bins = @(x) assert((length(x)>1) && isnumeric(x) && issorted(x), 'Bins must be monotonically increasing vector');
addParameter(p,'Bins',nan,valid_bins);
addParameter(p,'Axes',nan,@ishandle);
addParameter(p,'Smoothing','on', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Connect','on', @(x) any(validatestring(x,expectedFlags)));


% Parse inputs, save some to variables
parse(p,vects,places, varargin{:})
xspace = p.Results.XSpace;
bin_scale = p.Results.BinScale;
bins = p.Results.Bins;
ylim = p.Results.YLim;
colors = p.Results.Color;

% Create figure (if axes wasn't provided)
if ~ishandle(p.Results.Axes)
    viofig = figure('Position', [500, 1031, 800, 300], 'PaperPositionMode','auto');
    violin = axes('Parent',viofig);
else
    violin = p.Results.Axes;
end

% Get medians of all sets
medians = cellfun(@nanmedian,vects);

% Create figures; set XLim
if strcmp(p.Results.ShowBins,'on') % Diagnostic output: histogram overlaid with spline fit
    figure('Position',[500,357, 350 100*length(vects)]);
    ha = tight_subplot(length(vects),1);
end
if length(places)>1
    x_lim = [min(places)-(xspace*range(places)),max(places)+xspace*range(places)];
else
    x_lim = sort([places*0.9,places*1.1],'ascend');
end
tot_area = abs(diff(ylim(1:2)))*abs(diff(x_lim));

% Make bins
bin_scale = repmat(bin_scale,1,length(vects));
bin_scale = bin_scale(1:length(vects));

% Loop through sets, generate shapes
for i = 1:length(vects)
    if isnan(bins)
        bin_width = 2*iqr(all)*((numel(all)/length(vects))^(-1/3))/bin_scale(i);
        x = prctile(all,1):bin_width:prctile(all,99.5);
    else
        bin_width = bins(2)-bins(1);
        x = bins(:)';
    end
    
    
    % Generate histogram data
    y = hist(vects{i},x);
    
    % Cap histogram with zero values (keep spline from spiking @ end) and interpolate to get shape
    if y(1)==0
        pos = find(y>0,1,'first');
        y(1:(pos-1)) = [];
        x(1:(pos-1)) = [];
    end
    if y(end)==0
        pos = find(y>0,1,'last');
        y((pos+1):end) = [];
        x((pos+1):end) = [];
    end
    
    
    x = [min(x)-bin_width, x, max(x)+bin_width];
    y = [0 y/sum(y) 0];
    
    
    
    if strcmpi(p.Results.Smoothing,'on')
        xx = min(x):bin_width/10:max(x);
        yy = spline(x, y, xx);
        yy(yy<0) = 0;
    else
        xx = sort([x-bin_width/2.001, x+bin_width/2.001]);
        yy = repmat(y,2,1);
        yy = yy(:)';
    end
    
    % (optionally) show subplot of bins+spline fit
    if strcmp(p.Results.ShowBins,'on')
        hold(ha(i),'on')
        bar(ha(i),x,y,'FaceColor',[45 191 104]/255,'EdgeColor','none')
        set(ha(i),'XLim',ylim,'YLim',[0 .5])
        plot(ha(i),xx,yy,'LineWidth',2,'Color',[0 0 0])
        hold(ha(i),'off')
    end
    
    % Scale shape width so total area is consistent
    obj_width = p.Results.Area*tot_area/(sum(y)*diff(x(1:2))*2);

    % Make main violin plot
    if p.Results.LineWidth==0
        lnstyl = 'none';
        lnwid = 0.5;
    else
        lnstyl = '-';
        lnwid = p.Results.LineWidth;
    end
    
    hold(violin,'on')
    fill([places(i)+obj_width*yy,places(i)-obj_width*yy(end:-1:1)],[xx,xx(end:-1:1)],...
        colors{mod(i-1,length(colors))+1},'LineWidth',1,'Parent',violin,'LineStyle',lnstyl,'LineWidth',lnwid,...
        'EdgeColor',[0.4431 0.4510 0.4627])
    hold(violin,'off')
end
% Plot medians and set graph properties
hold(violin,'on')
if strcmpi(p.Results.Connect,'on')
    lnstyl = '-o';
else
    lnstyl = 'o';
end

plot(violin, places,medians,lnstyl,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
    'Color', [0 0 0],'LineWidth',1.2,'MarkerSize',7)
hold(violin,'off')
set(violin,'YLim',ylim,'XLim',x_lim);






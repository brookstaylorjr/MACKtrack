function violin = spaceviolin(vects, places, colors, graph_limits, show_bins, shape_area, bin_scale, xspace,axes_handle) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = spaceviolin(vects, places, colors, (graph_limits), (show_bins), (width), (bin_scale),(xspace),(axes_handle)) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SPACEVIOLIN makes a violin plot, spacing them according to a secondary vector (e.g. doses)
%
% 
% vects          1xN cell array of vectors (1-D array of object measurements)
% places         1xN array directing placement of each Violin
% colors         1x... cell vector specifying violin fill colors - cycles if length < N
% graph_limits   2 element vector with graph y-limits
% show_bins    (opt) create subplot showing histogram of each population, and spline-fit line
% shape_area   (opt) area-scaling for each violin (frac of total axis area) - default is 0.01
% bin_scale    (opt) scale num. of bins calculated by Freedman–Diaconis Rule - default=1
% xpace        (opt) how far to extend X axis limits graph past first/last place - % of range. Default = 0.1
% axes_handle  (opt) handle to axes where (new) violin figure will be plotted
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<9
    viofig = figure('Position', [500, 1031, 800, 300], 'PaperPositionMode','auto');
    violin = axes('Parent',viofig); % Violin-figure output
    if nargin<8
        xspace = 0.1;
        if nargin<7
            bin_scale = 1;
            if nargin<6
                shape_area = 0.01;
                if nargin<5
                    show_bins = 0;
                end
            end
        end
    end
else
    violin = axes_handle;
    viofig = get(violin,'Parent');
end

% Calculate bin sizes
all = cell2mat(vects(:));
IQR = iqr(all);
avg_size = length(all)/length(vects);
bin_width = 2*bin_scale*IQR/(avg_size^(1/3));
medians = zeros(size(vects));


% Make figures
if show_bins % Diagnostic output: histogram overlaid with spline fit
    figure('Position',[500,357, 350 100*length(vects)]);
    diagnos = tight_subplot(length(vects),1);
end

if length(places)>1
    x_lim = [min(places)-(xspace*range(places)),max(places)+xspace*range(places)];
else
    x_lim = [places*0.9,places*1.1];
end
tot_area = abs(diff(graph_limits(1:2)))*abs(diff(x_lim));


% Loop through sets
for i = 1:length(vects)
    % Generate histogram data
    x = min(vects{i}):bin_width:max(vects{i});
    h = hist(vects{i},x);
    medians(i) = nanmedian(vects{i});
    
    % Cap histogram with zero values (keep spline from spiking @ end) and spline-interpolate
    x = [min(x)-bin_width, x, max(x)+bin_width];
    y = [0 h/sum(h) 0];
    xx = min(x):bin_width/100:max(x);
    yy = spline(x, y, xx);
    yy(yy<0) = 0;
    
    % (optionally) show subplot of bins+spline fit
    if show_bins
        hold(diagnos(i),'on')
        bar(diagnos(i),x,y,'FaceColor',[45 191 104]/255,'EdgeColor','none')
        set(diagnos(i),'XLim',graph_limits,'YLim',[0 .5])
        plot(diagnos(i),xx,yy,'LineWidth',2,'Color',[0 0 0])
        hold(diagnos(i),'off')
    end
    
    % scale shape width so total area is consistent
    obj_width = shape_area*tot_area/(sum(yy)*bin_width/100*2);

    % Make main violin plot
    hold(violin,'on')
    fill([places(i)+obj_width*yy,places(i)-obj_width*yy(end:-1:1)],[xx,xx(end:-1:1)],...
        colors{mod(i-1,length(colors))+1},'LineWidth',2,'Parent',violin)
    hold(violin,'off')
end
% Plot medians and set graph properties
hold(violin,'on')
    plot(violin, places,medians,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
        'Color', [0 0 0],'LineWidth',2,'MarkerSize',8)
hold(violin,'off')
set(violin,'YLim',graph_limits,'XLim',x_lim);






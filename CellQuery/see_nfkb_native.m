function [graph, info, measure] = see_nfkb_native(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_NFKB_DIM is a data processing.and visualization script specialized to handle
% a nuclear-translocating species (it looks for NFkBdimNuclear and NFkBdimCytoplasm measurements).
%
% id             experiment ID (from Google Spreadsheet specigied in "loadID.m")
% show_graphs    boolean flag; specifies whether standard behavioral graphs will be shown
% diagnos        boolean flag; specifies whether optional diagnostic graphs will be shown
%
% 
% graph          primary output structure; must specify
%                   1) filtered/processed data (graph.var) 
%                   2) time vector for all images (graph.t) 
%                   3) XY convection adjustment (graph.shift) 
% info           secondary output structure; must specify
%                   1) Y limits for graphing (info.graph_limits)
%                   2) parameters from loadID.m (info.parameters) 
% measure         full output structure from loadID
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<3
    diagnos=0;
    if nargin<2
        show_graphs = 0;
    end
end

% Load data
[measure, info] = loadID(id);
info.Module = 'nfkbdimModule';

% Display parameters
max_shift = 2; % Max allowable frame shift in XY-specific correction
t_hrs = min([21,(size(measure.NFkBdimNuclear,2)-(1+2*max_shift))/12]); % Number of hours to display in graphs
info.graph_limits = [-0.5 8];
dendro = 0;
colors = setcolors;
robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr))));


% Filtering, part 1 cell fate and cytoplasmic intensity
droprows = [];
droprows = [droprows, sum(isnan(measure.NFkBdimCytoplasm(:,1:4)),2)>1]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:100)),2)>3]; % Long-lived cells
droprows = [droprows, sum(measure.NFkBdimCytoplasm(:,1:4)==0,2)>0]; % Very dim cells
%droprows = [droprows, info.CellData(:,end)]; % Non-edge cells

% NFkB normalization - subtract baseline for each cell; divide y background distribution width
nfkb = measure.NFkBdimNuclear(:,:);
nfkb_baseline = nanmin([prctile(nfkb(:,1:8),18.75,2),prctile(nfkb,10,2)],[],2);
nfkb = nfkb- repmat(nfkb_baseline,1,size(nfkb,2));
if diagnos
    figure,imagesc(nfkb,prctile(nfkb(:),[5,99])),colormap parula, colorbar
    title('All (baseline-subtracted) trajectories')
end

nfkb = nfkb/mean(info.parameters.adj_distr(2,:));


% Filtering, part 2: eliminate outlier cells (based on mean value)
nfkb_lvl = reshape(nfkb(max(droprows,[],2) == 0,:),[1 numel(nfkb(max(droprows,[],2) == 0,:))]);
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=3];
droprows =  [droprows, (nanmean(abs(nfkb-nanmean(nfkb_lvl)),2)./nanstd(nfkb_lvl))>=1.7];

% Filtering, part 3: nuclear stain intensity and starting NFkB value
keep = max(droprows,[],2) == 0;
start_lvl = prctile(nfkb(keep,1:8),18.75,2);
start_thresh = (nanmedian(start_lvl)+4*robuststd(start_lvl(:),2.5));
nuc_lvl = nanmedian(measure.MeanIntensityNuc(keep,1:31),2);
nuc_thresh = nanmedian(nuc_lvl)+2*robuststd(nuc_lvl(:),2);

droprows =  [droprows, prctile(nfkb(:,1:8),18.75,2) > start_thresh];
droprows =  [droprows, nanmedian(measure.MeanIntensityNuc(:,1:31),2) > nuc_thresh];

% Show some filter information
if diagnos
    filter_str = {'didn''t exist @ start', 'short-lived cells', 'NFkB<background',...
        'outliers [mean val >3*std]','extreme val [mean>1.7*std]', 'active @ start', 'high nuclear stain'};
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

    % Show cutoff for nuclear and starting level cutoffs 
    ranksmult(nfkb(keep,:),start_lvl);
    h = suptitle(['x = NFkB activation @ start. Threshold = ',num2str(start_thresh)]);
    set(h,'FontSize',14)
    ranksmult(nfkb(keep,:),nuc_lvl);
    h = suptitle(['x = Nuclear stain level. Threshold = ',num2str(nuc_thresh)]);
    set(h,'FontSize',14)


end

info.keep = max(droprows,[],2) == 0;
nfkb = nfkb(info.keep,:);

%% Initialize outputs, do final corrections

graph.t = 0:(1/info.parameters.FramesPerHour):t_hrs;
graph.celldata = info.CellData(info.keep,:);
graph.opt = maketicks(graph.t,info.graph_limits,0);
graph.opt.Name = 'NFkB Activation'; 
% Correct for XY positions that activate late
[graph.var, shift_xy] = alignTrajectories(nfkb(:,1:length(graph.t)+2*max_shift), graph.celldata, 60, max_shift);
[graph.order, graph.dendro.links] = hierarchial(graph.var(:,1:min([size(graph.var,2),150])),0);
if diagnos
    for i = 1:length(shift_xy)
        disp(['xy ',num2str(i),' shift : ',num2str(shift_xy(i))])
    end
end
graph.shift = shift_xy;




%% Graphing
if show_graphs
    % Colormap stack
    figs.a = figure('name','ColormapStack');
    set(figs.a,'Position', [500 700 1200 600])
    colormapStack(graph.var(graph.order,:),graph.celldata(graph.order,:), graph.opt);    
    
    % Hierarchial linkage
    if dendro
        graph.dendro.label = cell(size(graph.var,1),1);
        graph.dendro.label(:) = {''};
        figs.b = figure('name','EuclideanDist');
        graph.dendro.lines = dendrogram(graph.dendro.links,0,'orientation','left','labels',graph.dendro.label);
        set(figs.b,'Position', [555   743   300   535])
        set(get(graph.dendro.lines(1),'Parent'),'LooseInset',get(get(graph.dendro.lines(1),'Parent'),'TightInset')+ [0 0 0 0])
        set(graph.dendro.lines,'Color',[0.3 0.3 0.3])
    end

    % Line plot (mean+/-std)
    graph.line.top = nanmean(graph.var) + nanstd(graph.var);
    graph.line.bot = nanmean(graph.var) - nanstd(graph.var);
    figs.c = figure('name','MeanTrajectory');
    fill([graph.t,graph.t(end:-1:1)],[graph.line.top,graph.line.bot(end:-1:1)],colors.light_blue)
    hold on
    graph.line.main = plot(graph.t, nanmean(graph.var));
    set(graph.line.main,'Color',[.25 .25 .25],'LineWidth',2);
    hold off
    set(gca,'XTick',graph.opt.TimeTicks,'YTick',graph.opt.MeasurementTicks)
    set(gca,'YTickLabel',graph.opt.MeasurementTickLabels,'TickLength',[0.005 0.005])
    set(figs.c,'Position', [555   743   1145   300])
    axis([min(graph.opt.Times) max(graph.opt.Times) graph.opt.MeasurementBounds])

    % Small multiples line graph
    figs.d  = figure('name','smallmultiples');
    set(figs.d,'Position',[500, 350, 876, 1000]);
    graph.smult.order = randperm(size(graph.var,1),min([60 size(graph.var,1)]));
    graph.smult.h = tight_subplot(15,4);
    xpos = max(graph.opt.Times)-0.48*(max(graph.opt.Times)-min(graph.opt.Times));
    ypos =  max(graph.opt.MeasurementBounds) - 0.26*diff(graph.opt.MeasurementBounds);
    for i =1:length(graph.smult.order)
        plot(graph.smult.h(i),graph.t,graph.var(graph.smult.order(i),:),'Color',colors.grays{3}, 'LineWidth',2)
        set(graph.smult.h(i),'XLim',[min(graph.opt.Times) max(graph.opt.Times)],'YLim',graph.opt.MeasurementBounds)
        set(graph.smult.h(i),'XTickLabel',{[]}) 
        set(graph.smult.h(i),'YTickLabel',{[]}) 
        text(xpos,ypos,['XY ',num2str(graph.celldata(graph.smult.order(i),1)),...
            ', cell ',num2str(graph.celldata(graph.smult.order(i),2))],'Parent',graph.smult.h(i))
            
    end
end
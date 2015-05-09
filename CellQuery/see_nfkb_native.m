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
t_hrs = min([21,(size(measure.NFkBdimNuclear,2)-1)/12]); % Number of hours to display in graphs
max_shift = 0; % Max allowable frame shift in XY-specific correction
info.graph_limits = [-50 500];
dendro = 0;
colors = setcolors;

% Cytoplasmic trajectory -> stable, or is there some universal change (akin to RAW cells)?
if diagnos
    figure,ha = tight_subplot(2,1);
    plot(ha(1),smoothrows(measure.NFkBdimCytoplasm,7)')
    set(ha(1),'XTickLabel',{})
    plot(ha(2),nanmean(measure.NFkBdimCytoplasm))
    hold on
    plot(nanmedian(measure.NFkBdimCytoplasm),'Color',colors.blue)
    hold off
end

% Filter cells by fate
droprows = [];
droprows = [droprows, sum(isnan(measure.NFkBdimCytoplasm(:,1:4)),2)>1]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:100)),2)>3]; % Long-lived cells
droprows = [droprows, sum(measure.NFkBdimCytoplasm(:,1:4)==0,2)>0]; % Very dim cells
%droprows = [droprows, info.CellData(:,end)]; % Non-edge cells
info.keep = max(droprows,[],2) == 0;

% NFkB normalization
% Find baseline for each cell
nfkb = measure.NFkBdimNuclear(:,:);
nfkb_baseline = nanmin([prctile(nfkb(:,1:8),18.75,2),prctile(nfkb,10,2)],[],2);
nfkb = nfkb- repmat(nfkb_baseline,1,size(nfkb,2));
if diagnos
    figure,imagesc(nfkb,prctile(nfkb(:),[5,99])),colormap parula, colorbar
end
nfkb = nfkb(info.keep,:);


%% Filtering, part 2: eliminate outlier cells (based on mean value)
keep2 =  (nanmean(abs(nfkb-nanmean(nfkb(:))),2)./nanstd(nfkb(:)))<3;
info.keep(info.keep) =  keep2;
nfkb(~keep2,:) = [];

keep3 =  (nanmean(abs(nfkb-nanmean(nfkb(:))),2)./nanstd(nfkb(:)))<1.5;
info.keep(info.keep) =  keep3;
nfkb(~keep3,:) = [];


disp(['dropped ',num2str(sum(~keep2)),'+',num2str(+sum(~keep3)),' high-variance cells'])

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
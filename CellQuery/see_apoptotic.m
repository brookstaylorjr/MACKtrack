function [graph, info, measure] = see_apoptotic(id,show_graphs, diagnos)


%% Setup
if nargin<3
    diagnos=0;
    if nargin<2
        show_graphs = 0;
    end
end

%%

% Load data
[measure, info] = loadID(id);
info.Module = 'nfkbModule';

% Set display parameters
t_hrs = 12; % Number of hours to display in graphs
graph.t = 0:1/info.parameters.FramesPerHour:t_hrs;

max_shift = 0; % Max allowable frame shift in XY-specific correction
info.graph_limits = [-0.1 2];
dendro = 0;
colors = setcolors;

%% Filtering
bright = nanmin([prctile(measure.NFkBCytoplasm(:,1:8),18.75,2),nanmedian(measure.NFkBCytoplasm,2)],[],2);
bright_min = 0.9*info.parameters.nfkb_thresh;
bright_max = nanmean(bright(bright>bright_min)) + 1.8*nanstd(bright(bright>bright_min));

% Filter cells by fate and NFkB expression level
info.droprows = [];
info.droprows = [info.droprows, sum(isnan(measure.NFkBCytoplasm(:,1:4)),2)>1]; % Cells existing @ expt start
info.droprows = [info.droprows, sum(isnan(measure.NFkBNuclear(:,1:min([54,length(graph.t)]))),2)>3]; % Long-lived cells
%info.droprows = [info.droprows, info.CellData(:,end)]; % Non-edge cells
info.droprows = [info.droprows,(bright<bright_min)|(bright>bright_max)]; % Cells within brightness thresholds

% (graph of bright vs dim cells)
if diagnos
    brights_graph = measure.NFkBNuclear(max(info.droprows(:,1:end-2),[],2)==0,:);
    mean2 = bright(max(info.droprows(:,1:end-2),[],2)==0,:);
    [mean2,idx] = sort(mean2,'ascend');
    brights_graph = brights_graph(idx,:)./repmat(nanmean(brights_graph(idx,1:3),2),1,size(brights_graph,2));
    figure,ax = tight_subplot(3,1,0.05);
    imagesc(brights_graph(mean2<bright_min,:),'Parent',ax(1))
    set(ax(1),'CLim',[0.5 2.5],'XTick',[]), colormap gray
    title(ax(1),['Dim cells (starting nuclear level < ',num2str(round(bright_min*10)/10),')'])
    imagesc(brights_graph((mean2>bright_min)&(mean2<bright_max),:),'Parent',ax(2))
    set(ax(2),'CLim',[0.5 2.5],'XTick',[])
    imagesc(brights_graph(mean2>bright_max,:),'Parent',ax(3))
    set(ax(3),'CLim',[0.5 2.5],'XTick',[])
    title(ax(3),['Bright cells (starting nuclear level > ',num2str(bright_max),')'])
end

% Filter cells by nuclear/cytoplasmic expression levels
p1 = zeros(1,3);
p0 = zeros(1,3);
for i =1:3
    xvar1 = measure.NFkBCytoplasm(max(info.droprows,[],2) == 0,i);
    yvar1 = measure.NFkBNuclear(max(info.droprows,[],2) == 0,i);
    drop1 = isnan(xvar1)|isnan(yvar1);
    ply = polyfit(xvar1(~drop1),yvar1(~drop1),1);
    p1(i) = ply(1);
    p0(i) = ply(2);
end
nuc_adj = measure.NFkBNuclear-mean(p0);
ratio_early =  nanmedian(nuc_adj(:,1:3)./measure.NFkBCytoplasm(:,1:3),2);
info.droprows = [info.droprows,(ratio_early<(mean(p1)-0.2))...
    |(ratio_early>(mean(p1)+0.2))]; % Cells with physiological ratios

% Combine filters
info.keep = max(info.droprows,[],2) == 0;
%%
corr_nuc = measure.NFkBNuclear(info.keep,:);


corr_nuc = corr_nuc - repmat(prctile(corr_nuc(:,1:8),18.75,2),1,size(corr_nuc,2));
cyto_balanced = repmat(nanmedian(measure.NFkBCytoplasm(info.keep,1:5),2)-mean(info.parameters.img_distr(1,:))...
    ,1,size(corr_nuc,2));
corr_nuc = corr_nuc./cyto_balanced;

%% Filtering, part 2: do in 2 stages
keep2 =  (nanmean(abs(corr_nuc-nanmean(corr_nuc(:))),2)./nanstd(corr_nuc(:)))<3;
info.keep(info.keep) =  keep2;
corr_nuc(~keep2,:) = [];

keep3 =  (nanmean(abs(corr_nuc-nanmean(corr_nuc(:))),2)./nanstd(corr_nuc(:)))<1.8;
info.keep(info.keep) =  keep3;
corr_nuc(~keep3,:) = [];


disp(['dropped ',num2str(sum(~keep2)),'+',num2str(+sum(~keep3)),' high-variance cells'])

%% Initialize outputs, do final corrections

graph.celldata = info.CellData(info.keep,:);
graph.opt = maketicks(graph.t,info.graph_limits,0);
graph.opt.Name = 'NFkB Activation'; 

% Correct for XY positions that activate late
[graph.var, shift_xy] = alignTrajectories(corr_nuc(:,1:length(graph.t)+2*max_shift), graph.celldata, ...
    min([60,length(graph.t)]), max_shift);
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


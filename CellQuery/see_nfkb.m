function [graph, info, measure] = see_nfkb(id,show_graphs, diagnos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_NFKB is a data processing.and visualization script specialized to handle
% a nuclear-translocating species (it looks for NFkBNuclear and NFkBCytoplasm measurements).
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


%% Setup
if nargin<3
    diagnos=0;
    if nargin<2
        show_graphs = 0;
    end
end

% Load data
[measure, info] = loadID(id);
info.ImageExpr = info.parameters.nfkbModule.ImageExpr;

% Set display parameters
t_hrs = 12; % Number of hours to display in graphs
graph.t = 0:1/info.parameters.FramesPerHour:t_hrs;

max_shift = 0; % Max allowable frame shift in XY-specific correction
info.graph_limits = [-0.1 2.5];
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

%% Normalize dynamics for morphology changes and mixed expression levels 
% Identify and remove translocation-related artifacts in cytoplasm data
nfkb_cyto = measure.NFkBCytoplasm(info.keep,:);
nfkb_nuc = measure.NFkBNuclear(info.keep,:);
numcols = size(nfkb_nuc,2);
test1 = (nfkb_nuc - repmat(prctile(nfkb_nuc(:,1:8),18.75,2),1,numcols))./repmat(nanstd(nfkb_nuc,0,2),1,numcols);
test2 = (nfkb_cyto - repmat(nanmedian(nfkb_cyto(:,1:2),2),1,numcols))./repmat(nanstd(nfkb_cyto,0,2),1,numcols);
thresh = 2.5;
for i=1:size(nfkb_cyto,1)
    x1 = 1:size(nfkb_cyto,2);
    x1((test1(i,:)>thresh)&(test2(i,:)>thresh)) = [];
    y1 = nfkb_cyto(i,~((test1(i,:)>thresh)&(test2(i,:)>thresh)));
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];
    tmp = interp1(x1,y1,1:size(nfkb_cyto,2));
    tmp(isnan(nfkb_cyto(i,:))) = nan;
    nfkb_cyto(i,:) = tmp;
end

% Identify and remove data spikes from cytoplasm data
nfkb_cyto = nfkb_cyto./repmat(nanmean(nfkb_cyto(:,1:3),2),1,size(nfkb_cyto,2));
tmp = smoothrows(nfkb_cyto,5);
diffs = nfkb_cyto-tmp;
thresh = 5*nanstd(abs(diffs(:)));
deviat = abs(diffs)>thresh;
corr_cyto = nfkb_cyto;
for i=1:size(nfkb_cyto,1)
    x1 = find(~deviat(i,:));
    y1 = nfkb_cyto(i,~deviat(i,:));
    x1(isnan(y1)) = [];
    y1(isnan(y1)) = [];
    corr_cyto(i,:) = interp1(x1,y1,1:size(nfkb_cyto,2));
end
corr_cyto(isnan(nfkb_cyto)) = nan;

% Scale cytplasmic trajectories to balance nuclear baseline change.
smooth_cyto = smoothrows(corr_cyto,7);
corr_nuc = measure.NFkBNuclear(info.keep,:);
start_nuc = prctile(corr_nuc(:,1:8),18.75,2);
start_cyto = nanmedian(smooth_cyto(:,1:4),2);

win = 2;
std_nuc = stdfilt(corr_nuc,true(1,2*win+1));
std_nuc(imdilate(isnan(corr_nuc),true(1,2*win+1))) = inf;
omit_frames = ceil(0.8*info.parameters.FramesPerHour);

% Find each cell's late-time-point baseline and correct nuclear trajectory
for i =1:size(corr_nuc,1)
    min_v = min(std_nuc(i,omit_frames:end-win))*5;
    low_locs = find(std_nuc(i,:)<min_v);
    low_locs((low_locs<(omit_frames))|(low_locs>(size(std_nuc,2)-win))) = [];
    if numel(low_locs)>0
        % Find first point where cytoplasm and nuclear trajectories show same trend
        ind = 1;
        end_nuc = nanmean(corr_nuc(i,low_locs(ind)-win:low_locs(ind)+win));
        end_cyto = nanmean(smooth_cyto(i,low_locs(ind)-win:low_locs(ind)+win));
        while (((start_cyto(i)-end_cyto)/(start_nuc(i)-end_nuc)) < 0) && (ind<length(low_locs))
            ind = ind+1;
            end_nuc = nanmean(corr_nuc(i,low_locs(ind)-win:low_locs(ind)+win));
            end_cyto = nanmean(smooth_cyto(i,low_locs(ind)-win:low_locs(ind)+win));
        end
        % Find lowest baseline reached
        for j = (ind+1):length(low_locs)
            if ((start_cyto(i)-end_cyto)/(start_nuc(i)-end_nuc)) > 0
                nuc_a = nanmean(corr_nuc(i,low_locs(j)-win:low_locs(j)+win));
                cyto_a = nanmean(smooth_cyto(i,low_locs(j)-win:low_locs(j)+win));
                testval = nuc_a - abs(cyto_a-min([start_cyto(i),end_cyto]))...
                    /abs(start_cyto(i)-end_cyto)*abs(start_nuc(i)-end_nuc);
                if testval< end_nuc
                    end_nuc = nuc_a;
                    end_cyto = cyto_a;
                end
            end
        end
        % Apply correction (if nucleus and cytoplasm show same trend)
        if ((start_cyto(i)-end_cyto)/(start_nuc(i)-end_nuc)) > 0
            corr_row = (smooth_cyto(i,:) - min([start_cyto(i),end_cyto]))/abs(start_cyto(i)-end_cyto);
            corr_nuc(i,:) = (corr_nuc(i,:) - corr_row*abs(start_nuc(i)-end_nuc));
        end
        
    end
end

% Self-normalize corrected nuclear trajectories to starting cytoplasmic values
corr_nuc = corr_nuc - repmat(prctile(corr_nuc(:,1:8),18.75,2),1,size(corr_nuc,2));
cyto_balanced = repmat(nanmedian(measure.NFkBCytoplasm(info.keep,1:5),2)-mean(info.parameters.img_distr(1,:))...
    ,1,size(corr_nuc,2));
corr_nuc = corr_nuc./cyto_balanced;

% (graph of randomly selected uncorrected vs corrected trajectories)
if diagnos
    tmp_celldata = info.CellData(info.keep,:);
    nfkb_cyto = measure.NFkBCytoplasm(info.keep,:);
    nfkb_cyto = nfkb_cyto./repmat(nanmean(nfkb_cyto(:,1:3),2),1,size(nfkb_cyto,2));
    nfkb_nuc = measure.NFkBNuclear(info.keep,:);
    nfkb_nuc = nfkb_nuc./repmat(prctile(nfkb_nuc(:,1:8),18.75,2),1,size(nfkb_nuc,2));
    
    order1 = randperm(size(nfkb_cyto,1),6);
    legend1 = cell(size(order1));
    for i =1:length(order1)
        legend1{i} = ['xy',num2str(tmp_celldata(order1(i),1)),', cell ', num2str(tmp_celldata(order1(i),2))];
    end
    figure,ax = tight_subplot(4,1,0.08);
    plot(ax(1),nfkb_nuc(order1,:)','LineWidth',2),set(ax(1),'XLim',[1 size(nfkb_nuc,2)], 'YLim', [0.3 5])
    title(ax(1),'Raw nuclear NFkB (self-scaled)'), legend(ax(1),legend1)
    plot(ax(2),corr_nuc(order1,:)','LineWidth',2),set(ax(2),'XLim',[1 size(nfkb_nuc,2)], 'YLim', info.graph_limits)
    title(ax(2),'Corrected nuclear NFkB')
    plot(ax(3),nfkb_cyto(order1,:)','LineWidth',2),set(ax(3),'XLim',[1 size(nfkb_nuc,2)], 'YLim',[0.3 1.8])
    title(ax(3),'Raw cytoplasmic NFkB')
    plot(ax(4),smooth_cyto(order1,:)','LineWidth',2),set(ax(4),'XLim',[1 size(nfkb_nuc,2)], 'YLim', [0.3 1.8])
    title(ax(4),'Smoothed cytoplasmic NFkB')
end


%% Filtering, part 2: do in 2 stages
keep2 =  (nanmean(abs(corr_nuc-nanmean(corr_nuc(:))),2)./nanstd(corr_nuc(:)))<3;
info.keep(info.keep) =  keep2;
corr_nuc(~keep2,:) = [];

keep3 =  (nanmean(abs(corr_nuc-nanmean(corr_nuc(:))),2)./nanstd(corr_nuc(:)))<1.2;
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
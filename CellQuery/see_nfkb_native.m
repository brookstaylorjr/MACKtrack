function [graph, info, measure] = see_nfkb_native(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [graph, info, measure] = see_nfkb_native(id,graph_flag, verbose_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SEE_NFKB_NATIVE is a data processing.and visualization script specialized to handle
% a nuclear-translocating species (it looks for NFkBdimNuclear and NFkBdimCytoplasm measurements).
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'      'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'      'on' or 'off' - show verbose output
% 'MinLifetime'     final frame used to filter for long-lived cells (default = 100)
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
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);

% Optional parameters
expectedFlags = {'on','off'};
addParameter(p,'Display','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'Verbose','off', @(x) any(validatestring(x,expectedFlags)));
addParameter(p,'MinLifetime',100, @isnumeric);

% Parse parameters, assign to variables
parse(p,id, varargin{:})
if strcmpi(p.Results.Verbose,'on'); verbose_flag = 1; else verbose_flag = 0; end
if strcmpi(p.Results.Display,'on'); graph_flag = 1; else graph_flag = 0; end
MinLifetime = p.Results.MinLifetime;
%% Load data
[measure, info] = loadID(id);
info.Module = 'nfkbdimModule';

% Set display/filtering parameters
max_shift = 1; % Max allowable frame shift in XY-specific correction
start_thresh = 2; % Maximal allowable start level above baseline
info.graph_limits = [-0.25 8]; % Min/max used in graphing
dendro = 0;
colors = setcolors;


% Experiment-specific visualization settings/tweaks (load spreadsheet URL)
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')

% BT's experiments
if isnumeric(id)
    if  ~isempty(strfind(locations.spreadsheet,'10o_d9HN8dhw8bX4tbGxFBJ63ju7tODVImZWNrnewmwY'))
        % a) Heterozygous cell experiments
        if (id <= 270) || ismember(id,[370:379, 384:391, 395, 396])
            start_thresh = 1.5;
            info.graph_limits = [-0.25 5.5];
        end

        % b) bad light guide (4 experiments from same day)
        if ismember(id,315:318)
            info.graph_limits = [-0.2 4];
        end

        % c) 0.33ng TNF - delayed activation from slow mixing
        if id==290
            measure.NFkBdimNuclear = measure.NFkBdimNuclear(:,4:end);
            measure.NFkBdimCytoplasm = measure.NFkBdimNuclear(:,4:end);
            disp('Adjusted start point for this TNF expmt')
        end
        % d) 100uM CpG - delayed activation from slow mixing
        if id==283
            measure.NFkBdimNuclear = measure.NFkBdimNuclear(:,4:end);
            measure.NFkBdimCytoplasm = measure.NFkBdimNuclear(:,4:end);
            disp('Adjusted start point for this CpG expmt')
        end

        % e) prestimulated sets; don't filter pre-activated cells
        if ismember(id,[267:270, 323:324, 337:339,342:343, 356:360]); 
            start_thresh = 10;
        end

        % f) 09/23/2015 Missed an early timepoint; don't allow back-shifting
        if ismember(id,337:343) 
            max_shift = 0;
        end
        % g) IkBa KO set: poor nuclear staining (weak nuclei predominate): use eroded version of nuclear NFkB
        if ismember(id, 366:379)
            measure.NFkBdimNuclear = measure.NFkBdimNuclear_erode;
            disp('(Using eroded nuclei)')
        end
    % AA's experiments    
    elseif ~isempty(strfind(locations.spreadsheet,'1s51cOI7NvLIOEpEhZWjJKsPkCOE5Qz39m4tHa9nJ7ok'))
        % a) early experiments; heterozygous cells
        if id <= 60
            start_thresh = 1.5;
            info.graph_limits = [-0.25 6];
        else
            start_thresh = 1.8;
            info.graph_limits = [-0.25 5.5];
            info.baseline = 0.75;
        end
    end
end

%% Filtering
robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr))));

% Filtering, part 1 cell fate and cytoplasmic intensity
droprows = [];
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:4)),2)>2]; % Cells existing @ expt start
droprows = [droprows, sum(isnan(measure.NFkBdimNuclear(:,1:MinLifetime)),2)>3]; % Long-lived cells
droprows = [droprows, sum(measure.NFkBdimCytoplasm(:,1:4)==0,2)>0]; % Very dim cells
%droprows = [droprows, info.CellData(:,end)]; % Non-edge cells

% NFkB normalization - subtract baseline for each cell; divide y background distribution width
nfkb = measure.NFkBdimNuclear(:,:);
nfkb_smooth = nan(size(nfkb));
for i = 1:size(nfkb,1)
    nfkb_smooth(i,~isnan(nfkb(i,:))) = medfilt1(nfkb(i,~isnan(nfkb(i,:))),3);
end
% If default end frame is specified, use entire vector for baseline calculation. Otherwise use specified baseline.
if ismember('MinLifetime',p.UsingDefaults)
    nfkb_min = prctile(nfkb_smooth,2,2);
else
    nfkb_min = prctile(nfkb_smooth(:,1:MinLifetime),4,2);
end

nfkb_baseline = nanmin([nanmin(nfkb(:,1:4),[],2),nfkb_min],[],2);
nfkb = nfkb - repmat(nfkb_baseline,1,size(nfkb,2));
if verbose_flag
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
start_lvl = nanmin(nfkb(keep,1:3),[],2);
nuc_lvl = nanmedian(measure.MeanIntensityNuc(keep,1:31),2);
nuc_thresh = nanmedian(nuc_lvl)+2.5*robuststd(nuc_lvl(:),2);
area_thresh = 90;

droprows =  [droprows, prctile(nfkb(:,1:8),18.75,2) > start_thresh];
droprows =  [droprows, nanmedian(measure.MeanIntensityNuc(:,1:31),2) > nuc_thresh];

droprows =  [droprows, nanmedian(measure.Area,2) < area_thresh];

% Show some filter information
if verbose_flag
    filter_str = {'didn''t exist @ start', 'short-lived cells', 'NFkB<background',...
        'outliers [mean val >3*std]','extreme val [mean>1.7*std]', 'active @ start', 'high nuclear stain','low area'};
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
    
    ranksmult(nfkb(keep,:),nanmedian(measure.Area(keep,:),2))
    h = suptitle(['x = Median area. Threshold = ',num2str(area_thresh)]);
    set(h,'FontSize',14)

end

info.keep = max(droprows,[],2) == 0;
nfkb = nfkb(info.keep,:);

%% Initialize outputs, do final corrections
graph.celldata = info.CellData(info.keep,:);

% Correct for XY positions that activate late
[graph.var, shift_xy] = alignTrajectories(nfkb, graph.celldata, 60, max_shift);
if dendro
    [graph.order, graph.dendro.links] = hierarchial(graph.var(:,1:min([size(graph.var,2),150])),0);
else
    [~,graph.order] = sort(nansum(graph.var(:,1:min([size(graph.var,2),150])),2),'descend');
end
if verbose_flag
    for i = 1:length(shift_xy)
        disp(['xy ',num2str(i),' shift : ',num2str(shift_xy(i))])
    end
end
graph.shift = shift_xy;

graph.t = 0:(1/info.parameters.FramesPerHour):48;
graph.t = graph.t(1:min([length(graph.t),size(graph.var,2)]));
graph.opt = maketicks(graph.t,info.graph_limits,0);
graph.opt.Name = 'NFkB Activation'; 


%% Graphing
if graph_flag
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
        set(graph.smult.h(i),'XLim',[min(graph.opt.Times)-(range(graph.opt.Times)*0.02) max(graph.opt.Times)],...
            'YLim',graph.opt.MeasurementBounds,'XTickLabel',{},'YTickLabel',{})
        text(xpos,ypos,['XY ',num2str(graph.celldata(graph.smult.order(i),1)),...
            ', cell ',num2str(graph.celldata(graph.smult.order(i),2))],'Parent',graph.smult.h(i))
            
    end
end
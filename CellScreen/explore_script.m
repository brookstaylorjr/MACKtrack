% Sample script to examine some experimental results:


% First, load in the data
load('/Volumes/labdata/devon/Imaging experiments/20170201_TRELAP_adipogenesis2/AllData.mat')

% Load some other stuff, too
colors = setcolors; % define some pretty colors
colormaps = loadcolormaps; % define some pretty colormaps

%% - - - - - - - - SECTION 1: QUICKLY SUMMARIZING A FIXED CELL DATA EXPERIMENT - - - - - - - - - - - - - - - - 
summarizeMeasurement(AllData,'MeanNuc1') % Summary 1: look at distributions of a particular measurement across all conditions

conditions = getfields(AllData);
summarizeCondition(AllData.(conditions{end})) % Summary 2: can look at all measurements within a single condition -> make sure all images are ok (also shows all measurements made during analysis)

summarizeMeasurement2D(AllData,'MeanNuc1','MeanNuc2') % Summary 3: scatter plots of 2 variables of interest 


%% - - - - - - - - SECTION 2: REORGANIZE DATA to make it easier to pull out & compare selected conditions  - - - - - - - - 


[xdata_by_condition, xdata_by_well] = restructuredata(AllData,'MeanNuc2'); % 1st measurment we want to compare
[ydata_by_condition, ydata_by_well] = restructuredata(AllData,'MeanNuc1'); % A 2nd measurement we like to compare
condition_names = fieldnames(AllData); % Get names of experiments


subset = 9:16;
color_theme = colors.theme1(end:-1:1); % This is 4 colors: dark gray, dark blue, light blue, then red.
color_theme = repmat(color_theme,[1 ceil(length(subset)/length(color_theme))]); % Repeat colors if there are >4 subsets


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  EXAMPLE 2A: violin plots for selected experiments
% ~~~~~Modify this~~~~~~
x_data = xdata_by_condition;
% ~~~~~~~~~~~~~~~~~~~~~~~

% (Type 'help spaceviolin' to see the options you can set)
fig1 = figure('Position', positionfig(450, 390)); ax1 = axes('Parent',fig1);
spaceviolin(x_data(subset),1:length(subset),'Color',color_theme,'XSpace', 0.2,'Area',0.011,'YLim',[-500 6000],'Axes',ax1); 
ylabel('CEBP Expression'); 
set(gca,'XTick',1:4,'XTickLabel',{})
legend(condition_names(subset),'Location','northwest','Interpreter','none')

% VARIANT: compare multiple cell lines (log transform)
log_filter = @(vect) real(log(vect)); % (Make sure negative vals don't mess us up)
x_data = cellfun(log_filter,x_data,'UniformOutput',0);

subset2 = [5:8, 13:16, 21:24];
bins = linspace(5,10.5,40);
fig1 = figure('Position', positionfig(900, 200)); ax1 = axes('Parent',fig1);
spaceviolin(x_data(subset2),[1:4, 7:10, 13:16],'Color',color_theme,'XSpace', 0.2,'Area',0.005,'YLim',[4.75 11],...
    'Connect','off','Axes',ax1,'LineWidth',0,'Bins',bins); 
ylabel('CEBP Expression'); 
set(gca,'XTick',[],'XTickLabel',{})
legend(condition_names(subset2),'Location','northeast','Interpreter','none')


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2B: stacked/semi-transparent kernel density estimates (smoothed histograms)
ax1 = kdeoverlay(ydata_by_condition(subset),'Color',color_theme,'Alpha',0.25  ,'XLim',[-200, 3500],'LineWidth',2); 
xlabel('PPARG Expression'); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10,'Interpreter','none')

% Variation: do log transform on data 1st
log_filter = @(vect) real(log(vect)); % (Make sure negative vals don't mess us up)
x_data = cellfun(log_filter,xdata_by_condition,'UniformOutput',0);
[ax1, bw] = kdeoverlay(x_data(subset),'Color',color_theme,'Alpha',0.25, 'Xlim',[5 10.5],'LineWidth',2,'Bandwidth',0.12); 
xlabel('PPARG Expression (log)'); ylabel('Relative Frequency')



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2C: scatter plot of all conditions in subset (1 row of graphs)
% _____Modify these______
x_data = xdata_by_condition;
x_title = 'CEBPA expression';
y_data = ydata_by_condition;
y_title = 'PPARG expression';
% ________________________


figure('Position',positionfig(225*length(subset),235))
% Enforce consistent x- and y-scaling (using percentiles of the whole dataset)


% Automatically set graph limits
x_tmp = x_data{subset}; 
y_tmp = y_data{subset};
xlim = prctile(x_tmp(:),[0.01 99.999]);
ylim = prctile(y_tmp(:),[0.01 99.999]);
density_lim = [];

% Overrride with manually set graph and density limits
 xlim =[0 2000];
% ylim = [0 2000];
%density_lim = [0 5e-4]; % Sets the minimum/maximum density (i.e. color range) across all graphs


ha = tight_subplot(1,length(subset),[0.03 0.03],[0.2 0.1],[0.06 0.04],1);
for i = 1:length(subset)
    % Remove high outlier data
    x_tmp = x_data{subset(i)};
    y_tmp = y_data{subset(i)};
    outliers = (x_tmp>prctile(x_tmp,99.9)) | (y_tmp>prctile(y_tmp,99.9));
    
    
    [~,h] = dscatter2(x_data{subset(i)}(~outliers), y_data{subset(i)}(~outliers),'Parent',ha(i),'DensityLim',density_lim);
    alpha(h,0.05)
    set(ha(i),'XLim',xlim,'YLim',ylim,'XGrid','on','YGrid','on')
    text(mean(xlim),max(ylim),[num2str(subset(i)),') ', condition_names{subset(i)}],...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Parent',ha(i),'BackgroundColor','w',...
    'Interpreter','none')
    xlabel(ha(i), x_title)
    if i>1
        set(ha(i),'YTickLabel',{})
    else
        ylabel(ha(i),y_title)
    end
end
colormap(colormaps.viridis(end:-1:1,:))

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2D: Bar chart of (1) mean expression and (2) "differentiated cells" - percentage above some threshold.
% (Calculate these on a per well basis, including standard error)

% _____Modify these______
well_data = ydata_by_well;
threshold = 590;
name1 = 'Mean PPARG Expression';
name2 = 'Expressing cells';
xlabels = {'mCit 0dox', 'mCit D1dox', 'mCit D2dox', 'mCit D3dox', 'LAP 0dox', 'LAP D1dox', 'LAP D2dox', 'LAP D3dox'};
% ________________________


% Initialize structures to hold data
all_means = zeros(1, length(subset));
all_err =  zeros(1, length(subset));
pct_means = zeros(1, length(subset));
pct_err = zeros(1,length(subset));

% Fill in data
thresh_func = @(vect) sum(vect>threshold)/numel(vect);
for i = 1:length(subset)
    m = cellfun(@nanmean,well_data{subset(i)});
    all_means(i) = mean(m);
    all_err(i) = std(m)/sqrt(length(m));
    m = cellfun(thresh_func,well_data{subset(i)});
    pct_means(i) = mean(m);
    pct_err(i) = std(m)/sqrt(length(m));
end

% 1st figure: mean expression levels (w/ standard error, by well)
figure('Position',positionfig(400,200),'Name', 'Mean Values')
hold(gca,'on')
for i = 1:length(all_means)
    bar(i,all_means(i),'FaceColor',color_theme{i},'EdgeColor','none')
end
set(gca,'XTick',1:length(subset),'XTickLabel',xlabels,'XLim',[0 length(subset)+1])
h = terrorbar(1:length(all_means),all_means,all_err,0.4);
set(h,'Color',colors.grays{3},'LineWidth',2);
ylabel(name1,'FontWeight','bold')

% 2nd figure: mean percentage above threshold (w/ standard error, by well)
figure('Position',positionfig(400,200),'Name','Percentage above threshold')
hold(gca,'on')
for i = 1:length(pct_means)
    bar(i,pct_means(i),'FaceColor',color_theme{i},'EdgeColor','none')
end
set(gca,'XTick',1:length(subset),'XTickLabel',xlabels,'XLim',[0 length(subset)+1])
h = terrorbar(1:length(pct_means),pct_means,pct_err,0.4);
set(h,'Color',colors.grays{3},'LineWidth',2);
ylabel(name2,'FontWeight','bold')

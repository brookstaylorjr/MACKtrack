% Sample script to examine some experimental results:


% First, load in the data
load('/Volumes/labdata/devon/Imaging experiments/20170104_OP9TRE_doxcheck/Tracked/AllData.mat')

% Load some other stuff, too
colors = setcolors; % define some pretty colors
colormaps = loadcolormaps; % define some pretty colormaps

%% - - - - - - - - SECTION 1: QUICKLY SUMMARIZING A FIXED CELL DATA EXPERIMENT - - - - - - - - - - - - - - - - 
summarizeMeasurement(AllData,'MeanNuc1') % Summary 1: look at distributions of a particular measurement across all conditions

summarizeCondition(AllData.SD1LAP1dox) % Summary 2: can look at all measurements within a single condition -> make sure all images are ok. 

summarizeMeasurement2D(AllData,'Area','MeanNuc1') % Summary 3: scatter plots of 2 variables of interest 





%% - - - - - - - - SECTION2: REORGANIZE DATA to make it easier to pull out & compare selected conditions  - - - - - - - - 


[cebp_by_condition, cebp_by_well] = restructuredata(AllData,'MeanNuc1'); % 1st measurment we want to compare
[area_by_condition, area_by_well] = restructuredata(AllData,'Area'); % A 2nd measurement we like to compare
condition_names = fieldnames(AllData); % Get names of experiments


subset = 13:16; % Let's look at conditions #5-8,the dose response curve for SD1LAP0
color_theme = colors.theme1(end:-1:1); % This is 4 colors: dark gray, dark blue, light blue, then red.
color_theme = repmat(color_theme,[1 ceil(length(subset)/length(color_theme))]); % Repeat colors if there are >4 subsets


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  EXAMPLE 2A: violin plots for selected experiments
% (Type 'help spaceviolin' to see the options you can set)
fig1 = figure('Position', positionfig(450, 390)); ax1 = axes('Parent',fig1);
spaceviolin(cebp_by_condition(subset),1:length(subset),'Color',color_theme,'XSpace', 0.2,'Area',0.011,'YLim',[-500 6000],'Axes',ax1); 
ylabel('CEBP Expression'); 
set(gca,'XTick',1:4,'XTickLabel',{})
legend(condition_names(subset),'Location','northwest')


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2B: stacked/semi-transparent kernel density estimates (smoothed histograms)
ax1 = kdeoverlay(cebp_by_condition(subset),'Color',color_theme,'Alpha',0.25  ,'XLim',[-200, 3500],'LineWidth',2); 
xlabel('CEBP Expression'); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10)


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2C: scatter plot of all conditions in subset (1 row of graphs)
% _____Modify these______
x_data = cebp_by_condition;
x_title = 'CEPB expression';
y_data = area_by_condition;
y_title = 'Nuclear Area';
% ________________________


figure('Position',positionfig(800,200))
% Enforce consistent x- and y-scaling (using percentiles of the whole dataset)
xlim = prctile(cell2mat(x_data(:)),[.1 99]);
ylim = prctile(cell2mat(y_data(:)),[.1 99]);
ha = tight_subplot(1,length(subset),[0.03 0.03],[0.2 0.1],[0.06 0.04],1);
for i = 1:length(subset)
    dscatter(x_data{subset(i)}, y_data{subset(i)},'Parent',ha(i))
    set(ha(i),'XLim',xlim,'YLim',ylim,'XGrid','on','YGrid','on')
    text(mean(xlim),max(ylim),[num2str(i),') ', condition_names{i}],...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Parent',ha(i),'BackgroundColor','w')
    xlabel(ha(i), x_title)
    if i>1
        set(ha(i),'YTickLabel',{})
    else
        ylabel(ha(i),y_title)
    end
end
colormap(colormaps.byr)

%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2D: Bar chart of (1) mean expression and (2) "differentiated cells" - percentage above some threshold.
% (Calculate these on a per well basis, including standard error)

% _____Modify these______
x_data = cebp_by_well;
threshold = 595;
name1 = 'Mean CEBP Expression';
name2 = 'Expressing cells';
xlabels = {'0 dox', '0.5 dox', '1.0 dox', '1.5 dox'};
% ________________________


% Initialize structures to hold data
all_means = zeros(1, length(subset));
all_err =  zeros(1, length(subset));
pct_means = zeros(1, length(subset));
pct_err = zeros(1,length(subset));

% Fill in data
thresh_func = @(vect) sum(vect>threshold)/numel(vect);
for i = 1:length(subset)
    m = cellfun(@nanmean,x_data{subset(i)});
    all_means(i) = mean(m);
    all_err(i) = std(m)/sqrt(length(m));
    m = cellfun(thresh_func,x_data{subset(i)});
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

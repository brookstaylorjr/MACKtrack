% Sample script to examine HCS experimental results. (It's probably best to make a copy of this file, modify it as 
% required, and save it with your data).
% Brooks Taylor, 2/23/2017

% First, load in some data - chenage this file path to use a different "AllData" file.
%load('/Volumes/labdata/devon/Imaging experiments/20170206_OP9shCebpb_doxcheck/10h/AllData.mat')

% Load some other stuff, too
colors = setcolors; % define some pretty colors
colormaps = loadcolormaps; % define some pretty colormaps

% Summarize data contents
condition_names = fieldnames(AllData);
measurement_names = fieldnames(AllData.(condition_names{1}).Measurements);
disp('Conditions in this dataset:')
disp(condition_names)
disp('Measurements made (for cells in each condition):')
disp(measurement_names)
disp('- - - - - - - - ')

%% - - - - - - - - SECTION 1: QUICKLY SUMMARIZING A FIXED CELL DATA EXPERIMENT - - - - - - - - - - - - - - - - 
summarizeMeasurement(AllData,'MeanIntensity_cell2') % Summary 1: look at distributions of a particular measurement across all conditions


%%

summarizeCondition(AllData.(condition_names{1})) % Summary 2: look at all measurements within a SINGLE condition (divided by image)
summarizeMeasurement2D(AllData,'MeanNuc1','MeanNuc2') % Summary 3: scatter plots of 2 variables of interest 


%% - - - - - - - - SECTION 2: REORGANIZE DATA to make it easier to pull out & compare selected conditions  - - - - - - - - 

% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[xdata_by_condition, xdata_by_well] = restructuredata(AllData,'MeanNuc2'); % 1st measurment we want to compare
[ydata_by_condition, ydata_by_well] = restructuredata(AllData,'MeanNuc1'); % A 2nd measurement we want to compare
color_theme = colors.theme1(end:-1:1); % This is 4 colors: dark gray, dark blue, light blue, then red.
color_theme = repmat(color_theme,[1 ceil(length(condition_names)/length(color_theme))]); % Repeat colors if there are >4 subsets


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  EXAMPLE 2A: violin plots for selected experiments
% (Type 'help spaceviolin' to see the options you can set for this plot)


% ~~~~~~~~~~~   Parameters for violin plot 1 ~~~~~~~~~~~~~~~~~~
subset = 1:3; % Subset of conditions we want to plot
disp_data = ydata_by_condition(subset);
ylim = prctile(cell2mat(disp_data),[0.1 99.5]);
y_label1 = 'CEBPb expression';
area1 = 0.08; % Set width of each violin shape
padding = 0.5; % Set empty space at either end of violin plots
bin_scale = 0.6; % scale number of bins (DECREASE to make shape "smoother")
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(450, 390)); ax1 = axes('Parent',fig1);
spaceviolin(disp_data(subset),1:length(subset),'Color',color_theme,'XSpace', padding,'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1); 
ylabel(y_label1); 
set(gca,'XTick',1:length(subset),'XTickLabel',{})
legend(condition_names(subset),'Location','northwest','Interpreter','none')


% ~~~~~~~~~~~   Violin plot 2 (log transformed, larger/grouped subset of data) ~~~~~~~~~
subset2 = [5:8, 13:16]; % (Larger) set of subsets to plot
disp_data = ydata_by_condition(subset2);
groupings = [1:4, 7:10]; % Which groups of violin plots will be plotted together
area2 = 0.02; % Set width of each violin shape
y_label2 = 'CEBPb expression';
ylim2 = log(prctile(cell2mat(disp_data),[1 99.9]));
bins = linspace(ylim2(1),ylim2(2),25); % Manually set bins in each histogram here
padding2 = 0.3;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

log_filter = @(vect) real(log(vect)); % (Make sure negative vals don't mess us up)
disp_data = cellfun(log_filter,disp_data,'UniformOutput',0);
fig1 = figure('Position', positionfig(900, 200)); ax1 = axes('Parent',fig1);
spaceviolin(disp_data,groupings,'Color',color_theme,'XSpace', padding2,'Area',area2,'YLim',ylim2,...
    'Connect','off','Axes',ax1,'LineWidth',0,'Bins',bins); 
ylabel(y_label2); 
set(gca,'XTick',[],'XTickLabel',{})
legend(condition_names(subset2),'Location','northeast','Interpreter','none')


%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2B: overlaid semi-transparent kernel density estimates (smoothed histograms)

% ~~~~~~~~~~~~~~  Parameters for KDE overlay 1    ~~~~~~~~~~~~~
subset = 1:3;
disp_data = ydata_by_condition(subset);
alpha1 = 0.25; % Transparancy of each histogram
xlim = prctile(cell2mat(disp_data),[0.1 99.5]);
x_label1 = 'CEBPb expression';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ax1 = kdeoverlay(disp_data,'Color',color_theme,'Alpha',alpha1  ,'XLim',xlim,'LineWidth',2); 
xlabel(x_label1); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10,'Interpreter','none')


% ~~~~~~~~~~~~   Parameters for KDE overlay 2 (log-transformed) ~~~~~~~~~~~~~~~~~~
subset = 1:3;
disp_data = ydata_by_condition(subset);
alpha1 = 0.25; % Transparancy of each histogram
xlim2 = log(prctile(cell2mat(disp_data),[0.5 99.9]));
x_label1 = 'CEBPb expression';
bandwidth = 0.08; % Smoothing of KDE - can check MATLAB's guess by setting this to NaN and checking val of "bw", below

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

log_filter = @(vect) real(log(vect)); % (Make sure negative vals don't mess us up)
disp_data = cellfun(log_filter,disp_data,'UniformOutput',0);
[ax1, bw] = kdeoverlay(disp_data(subset),'Color',color_theme,'Alpha',0.25, 'Xlim',xlim2,'LineWidth',2); 
xlabel('PPARG Expression (log)'); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10,'Interpreter','none')



%% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EXAMPLE 2C: scatter plot of all conditions in subset (1 row of graphs)
% _____Modify these______
subset = 1:3;
x_data = xdata_by_condition(subset);
y_data = ydata_by_condition(subset);
x_title = 'CEBPA expression';
y_title = 'PPARG expression';
alpha1 = 0.25;

% Graph limits
xlim = prctile(cell2mat(x_data),[0.01 99]); % Can also use manual values here
ylim = prctile(cell2mat(y_data),[0.01 99]);
density_lim = []; % Leave empty if you want to automatically set.

% Overrride with manually set graph and density limits
% xlim =[0 2000];
% ylim = [0 2000];
density_lim = [0 5e-4]; % Sets the minimum/maximum density (i.e. color range) across all graphs
% ________________________


figure('Position',positionfig(225*length(subset),235))
ha = tight_subplot(1,length(subset),[0.03 0.03],[0.2 0.1],[0.06 0.04],1);
for i = 1:length(subset)
    % Remove high outlier data
    x_tmp = x_data{i};
    y_tmp = y_data{i};
    outliers = (x_tmp>prctile(x_tmp,99.9)) | (y_tmp>prctile(y_tmp,99.9));
    
    
    [~,h] = dscatter2(x_data{subset(i)}(~outliers), y_data{subset(i)}(~outliers),'Parent',ha(i),'DensityLim',density_lim);
    alpha(h,alpha1)
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
subset = 1:length(ydata_by_condition);
well_data = ydata_by_well(subset);
threshold = 2700;
name1 = 'Mean PPARg Expression';
name2 = 'Expressing cells';
xlabels = condition_names(subset);
% ________________________


% Initialize structures to hold data
all_means = zeros(1, length(subset));
all_err =  zeros(1, length(subset));
pct_means = zeros(1, length(subset));
pct_err = zeros(1,length(subset));

% Fill in data
thresh_func = @(vect) sum(vect>threshold)/numel(vect);
for i = 1:length(subset)
    m = cellfun(@nanmean,well_data{i});
    all_means(i) = mean(m);
    all_err(i) = std(m)/sqrt(length(m));
    m = cellfun(thresh_func,well_data{i});
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

% Sample script to examine HCS experimental results. (It's probably best to make a copy of this file, modify it as 
% required, and save it with your data).
% Brooks Taylor, 7/21/2017


% First, load in some data - change this file path to use a different "AllData" file.
load('/Volumes/labdata/brooks/Tracked/2016-10-24_physio-screen/Day 4/AllData.mat')
%load('Z:\brooks\Tracked\2016-10-24_physio-screen\Day 4\AllData.mat') % For Windows


% Summarize data contents
condition_names = fieldnames(AllData);
measurement_names = fieldnames(AllData.(condition_names{1}).Measurements);
disp('Conditions in this dataset:')
disp(condition_names)
disp('Measurements made (for cells in each condition):')
disp(measurement_names)
disp('- - - - - - - - ')

%% - - - - - - - - SECTION 1: QUICKLY SUMMARIZING A FIXED CELL DATA EXPERIMENT - - - - - - - - - - - - - - - - 
summarizeMeasurement(AllData,'MeanNuc1') % Summary 1: look at distributions of a particular measurement across all conditions

%%
summarizeMeasurement2D(AllData,'MeanNuc2','MeanNuc3') % Summary 2: scatter plots of 2 variables of interest

%%
summarizeCondition(AllData.OP9LD_DMI) % Summary 3: look at ALL measurements within a single condition - this can help identify low-quality images.


%% - - - - - - - - SECTION 2: FILTER AND REORGANIZE DATA  - - - - - - - - 

% Running these filters will ELIMINATE data points - don't uncomment these lines unless you want them!!

% (Optional) filter out poor-quality sites or wells from data
% AllData = dropImages(AllData,'_H05_s4');

% (Optional) filter out poor quality data (here, dim nuclei - likely false positives)
% thresh_func = @(dapi) (dapi < 550); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
% AllData = filterAllData(AllData,'MeanNuc1', thresh_func);

    
% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[xdata_by_condition, xdata_by_well] = restructuredata(AllData,'MeanNuc2'); % 1st measurment we want to compare/use
[ydata_by_condition, ydata_by_well] = restructuredata(AllData,'MeanNuc3'); % A 2nd measurement we want to compare/use



%% - - - - - - - - SECTION 3: EXAMPLES OF GRAPHING DATA  - - - - - - - - 
% After data is split by condition (or by well), you can easily pick a subset of conditions and compare them.


%% EXAMPLE 3A: a violin plot - type 'help violin' to see the full list of options you can set for this plot

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = xdata_by_condition; % Variable we want to plot
subset = 1:6; % Subset of conditions we want to plot
y_label = 'CEBPb expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 99]); % Y-limits of graph.
%ylim = [0 400]; % (Insead of automatically setting, can use manually-determined limits too)

area1 = 0.01; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 1; % Number of bins - decreaase to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(450, 390)); ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{})
legend(condition_names(subset),'Location','northwest','Interpreter','none')



%% EXAMPLE 3B: Smoothed histogram (KDE) overlay

% ~~~~~~~~~~~~~~  Parameters for KDE overlay  ~~~~~~~~~~~~~
disp_data = xdata_by_condition; % Variable we want to plot
subset = 1:10; % Subset of conditions we want to plot
x_label1 = 'CEBPb expression';

xlim = prctile(cell2mat(disp_data(subset)),[0.1 97]); % X-limits of graph
%xlim = [0 1000]; % (Insead of automatically setting, can use manually-determined limits too)

alpha1 = 0.3; % Transparancy of each histogram
bandwidth = nan; % Smoothing of KDE - can check MATLAB's guess by setting this to NaN and checking val of "bw", below
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[ax1, bw] = kdeoverlay(disp_data(subset),'Alpha',alpha1  ,'XLim',xlim,'LineWidth',2,'BandWidth',bandwidth); 
xlabel(x_label1); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10,'Interpreter','none')



%% EXAMPLE 3C: KDE overlay on log-transformed data

% ~~~~~~~~~~~~~~  Parameters for KDE overlay   ~~~~~~~~~~~~~
disp_data = xdata_by_condition; % Variable we want to plot
subset = 1:10; % Subset of conditions we want to plot
x_label1 = 'CEBPb expression';

xlim = log(prctile(cell2mat(disp_data(subset)),[0.1 99.5])); % X-limits of graph
%xlim = [0 1000]; % (Insead of automatically setting, can use manually-determined limits too)

alpha1 = 0.3; % Transparancy of each histogram
bandwidth = 0.11; % Smoothing of KDE - can check MATLAB's guess by setting this to NaN and checking val of "bw", below
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


log_filter = @(vect) real(log(vect)); % (Make sure negative vals don't mess us up)
disp_data = cellfun(log_filter,disp_data,'UniformOutput',0);
[ax1, bw] = kdeoverlay(disp_data(subset),'Alpha',alpha1, 'Xlim',xlim,'LineWidth',2); 
xlabel('PPARG Expression (log)'); ylabel('Relative Frequency')
legend(condition_names(subset),'Location','northeast','FontSize',10,'Interpreter','none')



%% EXAMPLE 2D: scatter plot of all conditions in subset (1 row of graphs)

% ~~~~~~~~~~~~~~  Parameters for scatter plot  ~~~~~~~~~~~~~
x_data = xdata_by_condition;
y_data = ydata_by_condition;
x_title = 'CEBPb expression';
y_title = 'PPARg expression';
subset = 1:8;

% Graph limits
xlim = prctile(cell2mat(x_data(subset)),[0.01 99]);
ylim = prctile(cell2mat(y_data(subset)),[0.01 99]);
% xlim =[0 2000]; % Override with manually set graph and density limits
% ylim = [0 2000];

alpha1 = 0.25; % Transparency of individual points
% density_lim = [0 5e-4]; % Sets the minimum/maximum density (i.e. color range) across all graphs
density_lim = []; % Leave empty if you want to automatically set.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


figure('Position',positionfig(225*length(subset),235))
ha = tight_subplot(1,length(subset),[0.03 0.03],[0.2 0.1],[0.06 0.04],1);
for i = 1:length(subset)
    % Remove high outlier data
    x_tmp = x_data{subset(i)};
    y_tmp = y_data{subset(i)};
    outliers = (x_tmp>prctile(x_tmp,99.9)) | (y_tmp>prctile(y_tmp,99.9))|isnan(x_tmp)|isnan(y_tmp);
    
    
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
colormap(cbrewer('seq','PuBu',256))


%% EXAMPLE 2E: Bar chart of (1) mean expression and (2) "differentiated cells" - percentage above some threshold.
% (Calculate these on a per well basis, including standard error)

% ~~~~~~~~~~~~~~  Parameters for bar plot  ~~~~~~~~~~~~~
well_data = xdata_by_well;
subset = 5:6; % Conditions you want to plot
name1 = 'Mean PPARg Expression';
name2 = 'Expressing cells';

threshold = 500; % Threshold to determine "percent expressing"
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Initialize some things.
colors = setcolors;
color_theme = colors.theme1(end:-1:1); % This is 4 colors: dark gray, dark blue, light blue, then red.
color_theme = repmat(color_theme,[1 ceil(length(condition_names)/length(color_theme))]); % Repeat colors if there are >4 subsets
xlabels = condition_names(subset);

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
set(gca,'XTick',1:length(subset),'XTickLabel',xlabels,'XLim',[0 length(subset)+1],'TickLabelInterpreter','None')
h = terrorbar(1:length(all_means),all_means,all_err,0.4);
set(h,'Color',colors.grays{3},'LineWidth',2);
ylabel(name1,'FontWeight','bold')

% 2nd figure: mean percentage above threshold (w/ standard error, by well)
figure('Position',positionfig(400,200),'Name','Percentage above threshold')
hold(gca,'on')
for i = 1:length(pct_means)
    bar(i,pct_means(i),'FaceColor',color_theme{i},'EdgeColor','none')
end
set(gca,'XTick',1:length(subset),'XTickLabel',xlabels,'XLim',[0 length(subset)+1],'TickLabelInterpreter','None')
h = terrorbar(1:length(pct_means),pct_means,pct_err,0.4);
set(h,'Color',colors.grays{3},'LineWidth',2);
ylabel(name2,'FontWeight','bold')

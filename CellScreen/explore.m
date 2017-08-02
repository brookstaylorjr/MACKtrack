% Sample script to examine HCS experimental results. (It's probably best to make a copy of this file, modify it as 
% required, and save it with your data).
% Brooks Taylor, 7/21/2017
% First, load in some data - change this file path to use a different "AllData" file.
load('/Volumes/labdata/kyle/Fixed_Cell_Data/170413 - MCF10A Hela Correlations MEK ERK GAPDH ENO MCM5 MCM7/Tracked3/MCF10A/AllData.mat')
% summarizeMeasurement(AllData,'IntegratedNuc4') % assess hists of Hoechst


%%
% Filter filter out poor quality data (here, dim nuclei - likely false positives)
thresh_func = @(nuc) (nuc < 2e4)|(nuc>7.5e4); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
AllData = filterAllData(AllData,'IntegratedNuc4', thresh_func);

%%
% Note: in this experiment, 'IntegratedCell' and 'MedianCell', and 'IntegratedCyto' and 'MedianCyto', are SWAPPED
% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[w2_by_condition, w2_by_well] = restructuredata(AllData,'IntegratedCell1'); % 1st measurment we want to compare/use
[w3_by_condition, w3_by_well] = restructuredata(AllData,'IntegratedCell2'); % A 2nd measurement we want to compare/use
[tot_by_condition, tot_by_well] = restructuredata(AllData,'IntegratedCell3'); % A 2nd measurement we want to compare/use


%%
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figure('Position',positionfig(1800,380))
ha = tight_subplot(2,9,[0.08 0.02],[0.08 0.03]);

for j = 1:numel(species_dist)
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,[99.9]), 64);
    hold(ha(j),'on')
    histogram(tmp,bins,'Parent',ha(j),'EdgeColor','none','FaceColor',colors.grays{3},'Normalization','pdf')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),200);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(j),x, g_pdf,'color',colors.navy)
    hold(ha(j),'off')
    
    
    set(ha(j),'XLim',[0 prctile(tmp,[99.9])],'XGrid', 'on','YGrid','on')
    
    text(mean(get(ha(j),'XLim')), max(get(ha(j),'YLim')),[species_list{j},' c.v. = ',...
        num2str(round(std(tmp) / mean(tmp)*100)/100) ],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(j),'FontSize',10)
end
disp('- - - - - ')
%% MCF10a - single species
select = [1 5 8 10 17 18];
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figs.hist_mcf10a = figure('Position',positionfig(1600,160));
ha = tight_subplot(1, length(select));

plot_idx = 1;
for j = select
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,99.8), 64);
    hold(ha(plot_idx),'on')
    histogram(tmp,bins,'Parent',ha(plot_idx),'EdgeColor','none','FaceColor',colors.grays{4},'Normalization','probability')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),120);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(plot_idx),x, g_pdf*(x(2)-x(1))*1.9,'color',colors.navy)
    hold(ha(plot_idx),'off')
    
    
    set(ha(plot_idx),'XLim',[min(bins) max(bins)],'XGrid', 'on','YGrid','on','YLim',[0 0.07])
    
    text(mean(get(ha(plot_idx),'XLim')), max(get(ha(plot_idx),'YLim')),...
        [species_list{j},' c.v. = ', num2str(round(std(tmp) / mean(tmp)*100)/100),...
        sprintf('\n'), ' \alpha =',num2str(a_b(1),'%.2g'),'  \beta =',num2str(a_b(2),'%.2g')],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(plot_idx),'FontSize',12)
 
    
    plot_idx = plot_idx+1;
    
end
disp('- - - - - ')





%% HELA CELLS

load('/Volumes/labdata/kyle/Fixed_Cell_Data/170413 - MCF10A Hela Correlations MEK ERK GAPDH ENO MCM5 MCM7/Tracked2/Hela/AllData.mat')
%%
summarizeMeasurement(AllData,'IntegratedNuc4') % assess hists of Hoechst


%%
% Filter filter out poor quality data (here, dim nuclei - likely false positives)
thresh_func = @(nuc) (nuc < 3e4)|(nuc>1.3e5); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
AllData = filterAllData(AllData,'IntegratedNuc4', thresh_func);

%%
% Note: in this experiment, 'IntegratedCell' and 'MedianCell', and 'IntegratedCyto' and 'MedianCyto', are SWAPPED
% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[w2_by_condition, w2_by_well] = restructuredata(AllData,'MedianCell1'); % 1st measurment we want to compare/use
[w3_by_condition, w3_by_well] = restructuredata(AllData,'MedianCell2'); % A 2nd measurement we want to compare/use
[tot_by_condition, tot_by_well] = restructuredata(AllData,'MedianCell3'); % A 2nd measurement we want to compare/use


%%
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figure('Position',positionfig(1800,380))
ha = tight_subplot(2,9,[0.08 0.02],[0.08 0.03]);

for j = 1:numel(species_dist)
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,[99.9]), 64);
    hold(ha(j),'on')
    histogram(tmp,bins,'Parent',ha(j),'EdgeColor','none','FaceColor',colors.grays{3},'Normalization','pdf')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),200);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(j),x, g_pdf,'color',colors.navy)
    hold(ha(j),'off')
    
    
    set(ha(j),'XLim',[0 prctile(tmp,[99.9])],'XGrid', 'on','YGrid','on')
    
    text(mean(get(ha(j),'XLim')), max(get(ha(j),'YLim')),[species_list{j},' c.v. = ',...
        num2str(round(std(tmp) / mean(tmp)*100)/100) ],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(j),'FontSize',10)
end
disp('- - - - - ')
%% MCF10a - single species
select = [1 5 8 10 17 18];
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figs.hist_hela = figure('Position',positionfig(1600,160));
ha = tight_subplot(1, length(select));

plot_idx = 1;
for j = select
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,99.8), 64);
    hold(ha(plot_idx),'on')
    histogram(tmp,bins,'Parent',ha(plot_idx),'EdgeColor','none','FaceColor',colors.lavender,'Normalization','probability')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),120);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(plot_idx),x, g_pdf*(x(2)-x(1))*1.9,'color',colors.navy)
    hold(ha(plot_idx),'off')
    
    
    set(ha(plot_idx),'XLim',[min(bins) max(bins)],'XGrid', 'on','YGrid','on','YLim',[0 0.07])
    
    text(mean(get(ha(plot_idx),'XLim')), max(get(ha(plot_idx),'YLim')),...
        [species_list{j},' c.v. = ', num2str(round(std(tmp) / mean(tmp)*100)/100),...
        sprintf('\n'), ' \alpha =',num2str(a_b(1),'%.2g'),'  \beta =',num2str(a_b(2),'%.2g')],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(plot_idx),'FontSize',12)
 
    
    plot_idx = plot_idx+1;
    
end
disp('- - - - - ')



%% _ - - - - - - - - - - - - - - - 


%%


colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figure('Position',positionfig(1800,560))
ha = tight_subplot(1,9,[0.08 0.02]);

for j = 1:length(species_dist)
    tmp_a = species_dist{1,j};
    tmp_b = species_dist{2,j};
    tmp2 = tot_by_condition{j};
    
    for i = 1; 
        cutoffs = [mean(tmp2)-3*std(tmp2), mean(tmp2)+3*std(tmp2)];
        tmp_a((tmp2<min(cutoffs)) | (tmp2>max(cutoffs))) = [];
        tmp_b((tmp2<min(cutoffs)) | (tmp2>max(cutoffs))) = [];

    end 
    disp(num2str(corr(tmp_a,tmp_b)))
    plot(ha(j),tmp_a,tmp_b,'.')
    
%     text(mean(get(ha(j),'XLim')), max(get(ha(j),'YLim')),[species_list{j},' c.v. = ',...
%         num2str(round(std(tmp) / mean(tmp)*100)/100) ],...
%         'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(j),'FontSize',10)
end

disp('- - - - - ')








%%




% Scan how CV changes as a function of outlier threshold (seems like 6-sigma is a pretty good cutoff)
for j = 1:9
    tmp = tot_by_condition{j};
    vect = linspace(3,12,100);
    hi_thresh = mean(tmp) + std(tmp)*vect;
    cvs = zeros(size(hi_thresh));
    len = zeros(size(hi_thresh));
    for i =1:length(hi_thresh); 
        cvs(i) = std(tmp(tmp<hi_thresh(i)))/ mean(tmp(tmp<hi_thresh(i)));
        len(i) = sum(tmp>hi_thresh(i));
    end
    figure,plotyy(vect,cvs,vect,len)
end



%%

for j = 1:9
    tmp = w2_by_condition{j};
    
    for i = 1:2
        cutoffs = [mean(tmp)-std(tmp)*5, mean(tmp)+5*std(tmp)];
        tmp(tmp<min(cutoffs)) = []; tmp(tmp>max(cutoffs)) = [];
    end
    
    
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    figure,histogram(tmp)
end
    
    
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


%%


% Sample script to examine HCS experimental results. (It's probably best to make a copy of this file, modify it as 
% required, and save it with your data).
% Brooks Taylor, 7/21/2017


% First, load in some data - change this file path to use a different "AllData" file.
load('/Volumes/labdata/kyle/Fixed_Cell_Data/170413 - MCF10A Hela Correlations MEK ERK GAPDH ENO MCM5 MCM7/Tracked2/MCF10A/AllData.mat')
%%
summarizeMeasurement(AllData,'IntegratedNuc4') % assess hists of Hoechst


%%
% Filter filter out poor quality data (here, dim nuclei - likely false positives)
thresh_func = @(nuc) (nuc < 2e4)|(nuc>7.5e4); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
AllData = filterAllData(AllData,'IntegratedNuc4', thresh_func);

%%
% Note: in this experiment, 'IntegratedCell' and 'MedianCell', and 'IntegratedCyto' and 'MedianCyto', are SWAPPED
% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[w2_by_condition, w2_by_well] = restructuredata(AllData,'MedianCell1'); % 1st measurment we want to compare/use
[w3_by_condition, w3_by_well] = restructuredata(AllData,'MedianCell2'); % A 2nd measurement we want to compare/use
[tot_by_condition, tot_by_well] = restructuredata(AllData,'MedianCell3'); % A 2nd measurement we want to compare/use


%%
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figure('Position',positionfig(1800,380))
ha = tight_subplot(2,9,[0.08 0.02],[0.08 0.03]);

for j = 1:numel(species_dist)
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,[99.9]), 64);
    hold(ha(j),'on')
    histogram(tmp,bins,'Parent',ha(j),'EdgeColor','none','FaceColor',colors.grays{3},'Normalization','pdf')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),200);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(j),x, g_pdf,'color',colors.navy)
    hold(ha(j),'off')
    
    
    set(ha(j),'XLim',[0 prctile(tmp,[99.9])],'XGrid', 'on','YGrid','on')
    
    text(mean(get(ha(j),'XLim')), max(get(ha(j),'YLim')),[species_list{j},' c.v. = ',...
        num2str(round(std(tmp) / mean(tmp)*100)/100) ],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(j),'FontSize',10)
end
disp('- - - - - ')
%% MCF10a - single species
select = [1 5 8 10 17 18];
colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figs.hist_mcf10a = figure('Position',positionfig(1600,160));
ha = tight_subplot(1, length(select));

plot_idx = 1;
for j = select
    tmp = species_dist{j};
    for i = 1:3;
        cutoffs = [mean(tmp)-6*std(tmp), mean(tmp)+6 *std(tmp)];
        tmp((tmp<min(cutoffs)) | (tmp>max(cutoffs))) = [];
    end 
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    bins = linspace(0, prctile(tmp,99.8), 64);
    hold(ha(plot_idx),'on')
    histogram(tmp,bins,'Parent',ha(plot_idx),'EdgeColor','none','FaceColor',colors.grays{4},'Normalization','probability')
    
    % Fit and plot a gamma distribution
    a_b = gamfit(tmp);
    x = linspace(min(bins), max(bins),120);
    g_pdf = gampdf(x,a_b(1),a_b(2));
    plot(ha(plot_idx),x, g_pdf*(x(2)-x(1))*1.9,'color',colors.navy)
    hold(ha(plot_idx),'off')
    
    
    set(ha(plot_idx),'XLim',[min(bins) max(bins)],'XGrid', 'on','YGrid','on','YLim',[0 0.07])
    
    text(mean(get(ha(plot_idx),'XLim')), max(get(ha(plot_idx),'YLim')),...
        [species_list{j},' c.v. = ', num2str(round(std(tmp) / mean(tmp)*100)/100),...
        sprintf('\n'), ' \alpha =',num2str(a_b(1),'%.2g'),'  \beta =',num2str(a_b(2),'%.2g')],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(plot_idx),'FontSize',12)
 
    
    plot_idx = plot_idx+1;
    
end
disp('- - - - - ')




%%


colors = setcolors;
species_list = {'MCM7', 'GAPDH','GAPDH','MEK1', 'MEK1/2', 'GAPDH','GAPDH','GAPDH','GAPDH',;
                'MCM5', 'MCM7', 'MCM5', 'ERK', 'ERK',    'MEK1',  'ERK', 'Eno1', 'MEK1/2'}';
species_dist = [w2_by_condition, w3_by_condition]';
tot_dist = [tot_by_condition, tot_by_condition]';

% [species_list,order] = sort(species_list);
% species_dist = species_dist(order);

figure('Position',positionfig(1800,560))
ha = tight_subplot(1,9,[0.08 0.02]);

for j = 1:length(species_dist)
    tmp_a = species_dist{1,j};
    tmp_b = species_dist{2,j};
    tmp2 = tot_by_condition{j};
    
    for i = 1; 
        cutoffs = [mean(tmp2)-3*std(tmp2), mean(tmp2)+3*std(tmp2)];
        tmp_a((tmp2<min(cutoffs)) | (tmp2>max(cutoffs))) = [];
        tmp_b((tmp2<min(cutoffs)) | (tmp2>max(cutoffs))) = [];

    end 
    disp(num2str(corr(tmp_a,tmp_b)))
    plot(ha(j),tmp_a,tmp_b,'.')
    
%     text(mean(get(ha(j),'XLim')), max(get(ha(j),'YLim')),[species_list{j},' c.v. = ',...
%         num2str(round(std(tmp) / mean(tmp)*100)/100) ],...
%         'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(j),'FontSize',10)
end

disp('- - - - - ')








%%




% Scan how CV changes as a function of outlier threshold (seems like 6-sigma is a pretty good cutoff)
for j = 1:9
    tmp = tot_by_condition{j};
    vect = linspace(3,12,100);
    hi_thresh = mean(tmp) + std(tmp)*vect;
    cvs = zeros(size(hi_thresh));
    len = zeros(size(hi_thresh));
    for i =1:length(hi_thresh); 
        cvs(i) = std(tmp(tmp<hi_thresh(i)))/ mean(tmp(tmp<hi_thresh(i)));
        len(i) = sum(tmp>hi_thresh(i));
    end
    figure,plotyy(vect,cvs,vect,len)
end



%%

for j = 1:9
    tmp = w2_by_condition{j};
    
    for i = 1:2
        cutoffs = [mean(tmp)-std(tmp)*5, mean(tmp)+5*std(tmp)];
        tmp(tmp<min(cutoffs)) = []; tmp(tmp>max(cutoffs)) = [];
    end
    
    
    disp(['cond. ' ,num2str(j),' c.v.: ',num2str(std(tmp) / mean(tmp))])
    figure,histogram(tmp)
end
    
    
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


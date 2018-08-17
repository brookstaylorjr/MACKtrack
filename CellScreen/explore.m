% Testing with 2 dCas9-based constructs
% - a KRAB-dCas9-DHFR
% - a Cas9-P2A-VP64-PP7-DHFR


% I did a single dose response (0, 200, 2000 uM TMP) , measuring protein expression 18 hrs after induction. sgRNA 
% constructs targeted Cebpb or Pparg2 (3 variants for each gene)

% First, load in some data - change this file path to use a different "AllData" file.
load('/Volumes/labdata/zahra/2018 Fixed data/20180112-ZB-LIN/Tracked/Lin-First/AllData.mat')
%load('Z:\brooks\Tracked\2016-10-24_physio-screen\Day 4\AllData.mat') % For Windows

%%
% Summarize data contents
condition_names = fieldnames(AllData);
measurement_names = fieldnames(AllData.(condition_names{1}).Measurements);
disp('Conditions in this dataset:')
disp(condition_names)
disp('Measurements made (for cells in each condition):')
disp(measurement_names)
disp('- - - - - - - - ')



%% - - - - - - - - FILTER AND REORGANIZE DATA  - - - - - - - - 

% Running these filters will ELIMINATE data points - don't uncomment these lines unless you want them!!

% (Optional) filter out poor-quality sites or wells from data
% AllData = dropImages(AllData,'_H05_s4');

% (Optional) filter out poor quality data (here, dim nuclei - likely false positives)
thresh_func = @(dapi) (dapi >4.5e6)|(dapi<1e6); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
AllData = filterAllData(AllData,'IntegratedNuc1', thresh_func);

% Some issues with poor fluorescence on far-away cells - filter them out.
threshold_radius = @(x,y) (sqrt((1080 - x).^2 + (1080 - y).^2))>1100;
AllData = filterAllData2(AllData,'CentroidX','CentroidY', threshold_radius);

    
% Peel off/restructure 2 data fields of interest (e.g. PPARg and CEBPb intensity)
[pparg_by_condition, pparg_by_well] = restructuredata(AllData,'MeanNuc2'); % 1st measurment we want to compare/use
[cebpb_by_condition, cebpb_by_well] = restructuredata(AllData,'MeanNuc3'); % A 2nd measurement we want to compare/use

summarizeMeasurement(AllData,'IntegratedNuc1') % Summary 1: look at distributions of a particular measurement across all conditions


%% 1) C/EBPb expression in "activator" constructs

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = cebpb_by_condition; % Variable we want to plot
subset = 26:38; % Subset of conditions we want to plot
y_label = 'C/EBPb expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 97]); % Y-limits of graph.

area1 = 0.012; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 1; % Number of bins - decrease to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(890, 240)); 

ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1,'Connect','off'); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{},'XLim',[0 length(subset)+5])
legend(condition_names(subset),'Location','east','Interpreter','none')


%% 1) C/EBPb expression in "activator" constructs

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = cebpb_by_condition; % Variable we want to plot
subset = [26 20 27 33 21 28 34 22 29 35 23 30 36 24 31 37 25 32 38]; % Subset of conditions we want to plot
y_label = 'C/EBPb expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 97.5]); % Y-limits of graph.

area1 = 0.012; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 0.95; % Number of bins - decreaase to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(890, 240)); 

ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1,'Connect','off'); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{},'XLim',[0 length(subset)+5])
legend(condition_names(subset),'Location','east','Interpreter','none')


%% 2) C/EBPb expression in "inhibitor" constructs

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = cebpb_by_condition; % Variable we want to plot
subset = [7 1 8 14 2 9 15 3 10 16 4 11 17 5 12 18 6 13 19]; % Subset of conditions we want to plot
y_label = 'C/EBPb expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 97.5]); % Y-limits of graph.

area1 = 0.012; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 0.95; % Number of bins - decreaase to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(890, 240)); 

ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1,'Connect','off'); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{},'XLim',[0 length(subset)+5])
legend(condition_names(subset),'Location','east','Interpreter','none')

%% 3) PPARg expression in "activator" constructs

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = pparg_by_condition; % Variable we want to plot
subset = [26 20 27 33 21 28 34 22 29 35 23 30 36 24 31 37 25 32 38]; % Subset of conditions we want to plot
y_label = 'PPARg expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 97.5]); % Y-limits of graph.

area1 = 0.012; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 0.95; % Number of bins - decreaase to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(890, 240)); 

ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1,'Connect','off'); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{},'XLim',[0 length(subset)+5])
legend(condition_names(subset),'Location','east','Interpreter','none')

%% 4) PPARg expression in "inhibitor" constructs

% ~~~~~~~~~~~   Basic parameters for a violin plot  ~~~~~~~~~~~~~~~~~~
disp_data = pparg_by_condition; % Variable we want to plot
subset = [7 1 8 14 2 9 15 3 10 16 4 11 17 5 12 18 6 13 19]; % Subset of conditions we want to plot
y_label = 'PPARg expression'; % Y-label of graph

ylim = prctile(cell2mat(disp_data(subset)),[0.01 97.5]); % Y-limits of graph.

area1 = 0.012; % Sets width of each violin shape - increase to make shape FATTER.
bin_scale = 0.95; % Number of bins - decreaase to make shape SMOOTHER.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig1 = figure('Position', positionfig(890, 240)); 

ax1 = axes('Parent',fig1);
violin(disp_data(subset),1:length(subset), 'Area',area1,'YLim',ylim,....
    'BinScale',bin_scale, 'Axes',ax1,'Connect','off'); 
ylabel(y_label); 
set(gca,'XTick',1:length(subset),'XTickLabel',{},'XLim',[0 length(subset)+5])
legend(condition_names(subset),'Location','east','Interpreter','none')


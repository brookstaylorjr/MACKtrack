% Sample script to examine HCS experimental results. (It's probably best to make a copy of this file, modify it as 
% required, and save it with your data).
% Brooks Taylor, 7/21/2017


% First, load in some data - change this file path to use a different "AllData" file.
load('/Volumes/labdata/brooks/Tracked/2017-11-10_FBS-test/Day 2/AllData.mat')
AllDays{1}= AllData;
load('/Volumes/labdata/brooks/Tracked/2017-11-10_FBS-test/Day 4/AllData.mat')
AllDays{2}= AllData;


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
summarizeMeasurement(AllDays{2},'IntegratedNuc1') % Summary 1: look at distributions of a particular measurement across all conditions

%%
summarizeMeasurement2D(AllData,'MeanNuc2','MeanNuc3') % Summary 2: scatter plots of 2 variables of interest

%%
summarizeCondition(AllData.OP9LD_DMI) % Summary 3: look at ALL measurements within a single condition - this can help identify low-quality images.


%% - - - - - - - - SECTION 2: FILTER AND REORGANIZE DATA  - - - - - - - - 

% Running these filters will ELIMINATE data points - don't uncomment these lines unless you want them!!

% (Optional) filter out poor-quality sites or wells from data
% AllData = dropImages(AllData,'_H05_s4');

% (Optional) filter out poor quality data (here, dim nuclei - likely false positives)
thresh_func = @(dapi) (dapi < 1.6e6) | (dapi > 5.5e6); % (you can pick a value by inspecting histograms from "summarizeMeasurement")
for i = 1:length(AllDays)
   AllDays{i} = filterAllData(AllDays{i},'IntegratedNuc1', thresh_func);
end
%%


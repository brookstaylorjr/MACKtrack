% Single layout scans (define measurements, directory, and parameters, then call buildFromLayout)
% Build measurment sets from layout files

%% Prat's 4-plate degradation/synthesis inhibition experiment - April 8th, 2015
% Specify staring directory - we'll find ALL layouts in any subdirectories under this one
start_dir = 'Z:\prat\20160428 Experiment';
% Specify channels and measurements
nuclear_channel = 'w1';
measurement_channels = {'w2','w1'};
measurement_type = {@mean,@mean};
% Load parameters
load('OP9_default.mat')
% Find all plate layouts in parent directory; segment and measure cells
[ScreenData, layout_dir, image_dir] = ...
    buildFromLayout(start_dir, nuclear_channel, measurement_channels, measurement_type, parameters);
% Save output file(s)
save([start_dir,filesep,'ScreenData_',date,'.mat'],'ScreenData')
disp('... saved file')



%%

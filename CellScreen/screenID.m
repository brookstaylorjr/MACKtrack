function [] = screenID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = screenID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCREEN queues and tracks multiple experimental sets from the "Scope Runs" Google Doc -  
% choose sets by ID number. Similar to runID in functionality, but uses associated 
% layout.xlsx files (saved with images) to identify, group, segment, and measure cells.
%
% INPUT:
% varargin     ID# of sets to track (1st column of spreadsheet) - single values or vector
%
% Example Usage: sceenID(1,4:16,27) will sequentially analyze sets 1, 4 to 16, and 27.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Specify save directory and load spreadsheet URL
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end
[~,status] = urlread(locations.spreadsheet);
if ~status
    error(['Spreadsheet URL load unsuccessful. Please load, update, and re-save "locations.mat"'])
end

% Read in content from the "Scope Runs - Tracking Sets" spreadsheet
data = readScopeRuns(locations.spreadsheet, cell2mat(varargin));


% Display the sets we're running
for idx = 1:numel(data.save_folder)
    disp(['- ',data.save_folder{idx}])
end




% Cycle/measure sets
for idx = 1:numel(data.save_folder)
    
    % PARAMETERS
    load([home_folder,'Parameters',filesep,data.parameter_files{idx}])
    
    [ScreenData, layout_dir, image_dir] = ...
    buildFromLayout(start_dir, nuclear_channel, measurement_channels, measurement_type, parameters);
    
    
end
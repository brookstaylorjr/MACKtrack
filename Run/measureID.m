function [] = measureID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MEASUREID queues and measures multiple experimental sets from the "Scope Runs" Google Doc 
% Choose sets by ID number. This function is identical to runID, but assumes cells have been
% tracked already.
%
% ALL DATA IS SAVED IN '/home/brooks/Data'
%
% varargin     ID# of sets to track (1st column of spreadsheet)
%
% Example Usage: runById(1,4:16,27) will track, sequentially, 15 experimental sets using
% the parameters specified in the Google Doc.
% 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Specify save directory and load spreadsheet URL
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
[~,status] = urlread(locations.spreadsheet);
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end
if ~status
    error(['Spreadsheet URL load unsuccessful. Please load, update, and re-save "locations.mat"'])
end


% Read in content from the "Scope Runs - Tracking Sets" spreadsheet
data = readScopeRuns(locations.spreadsheet, cell2mat(varargin));


% Display the sets we're running
for idx = 1:numel(data.dates)
    disp(['- ',data.dates{idx},'_',data.names{idx}])
end

% Cycle/measure sets
for idx = 1:numel(data.dates)
    
    % PARAMETERS
    load([home_folder,'Parameters',data.parameter_files{idx}])
    parameters.ImagePath = data.image_paths{idx};
    parameters.TimeRange = eval(data.time_ranges{idx});
    parameters.XYRange = eval(data.xy_ranges{idx});
    parameters.SaveDirectory = [data.save_dir{idx},filesep,data.dates{idx},'_',data.names{idx}];
    
    
    % MEASUREMENT
    disp('Measuring...')
    try
        UCSDcellMeasure(parameters);      
    catch ME
        disp(['Error in measurement:' , ME.message])
        for err = 1:length(ME.stack)
            disp(['-> ', ME.stack(err).name,', line ', num2str(ME.stack(err).line)])
        end
        error('killing measurement')
    end
end
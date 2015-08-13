function [] = runID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = runID(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% RUNID queues and tracks multiple experimental sets from the "Scope Runs" Google Doc -  
% choose sets by ID number.
%
% INPUT:
% varargin     ID# of sets to track (1st column of spreadsheet) - single values or vector
%
% Example Usage: runID(1,4:16,27) will sequentially track sets 1, 4 to 16, and 27.
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
    parameters.ImagePath = data.image_paths{idx};
    parameters.TimeRange = eval(data.time_ranges{idx});
    parameters.XYRange = eval(data.xy_ranges{idx});
    parameters.SaveDirectory = [data.save_dir{idx},filesep,data.save_folder{idx}];
    clear p;
    eval(data.modify{idx});
    if exist('p','var'); parameters = combinestructures(p,parameters); end;

    mkdir([locations.data,filesep,parameters.SaveDirectory])
    % TRACKING
    parfor i = 1:length(parameters.XYRange)
        xyPos = parameters.XYRange(i);
        try
            if strcmp(parameters.ImageType,'None')
                trackPrimary(parameters,xyPos)
            else
                trackLoop(parameters,xyPos) % DIC or phase
            end
        catch ME
            disp(['Error in tracking position ', num2str(xyPos),':' , ME.message])
            for err = 1:length(ME.stack)
                disp(['-> ', ME.stack(err).name,', line ', num2str(ME.stack(err).line)])
            end
            error('killing tracking')
        end
    end
    
    % MEASUREMENT
    disp(['Measuring ', parameters.SaveDirectory,'...'])
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
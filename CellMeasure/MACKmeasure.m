function [] = MACKmeasure(parameters,parallel_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = MACKmeasure(parameters,parallel_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MACKMEASURE calculates morphological features and optional
% additional-channel data from a MACKtrack output set.
%
% Takes in a MACKtrack output directory, loops through all xy positions,
% creating individual measurement sets (using the primarymeasure
% subfunction) for each.
% 
% At the end, it concatenates measurement sets into AllMeasurements 
% structure, which is saved in the output directory
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<2
    parallel_flag = parameters.Parallel;
end

% Make directory names, list contents, initialize variables 
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
parameters.locations = locations;
parameters.TotalImages = length(parameters.TimeRange);
orig_params = parameters; % Save params before they're flatfield-modified

% Convert any parameter flatfield images to functions; add background image
if isfield(parameters,'Flatfield')
    parameters.Flatfield = processFlatfields(parameters.Flatfield);
end


% Outer loop: Cycle xy folders (in each condition)
if parallel_flag
    parfor i = 1:length(parameters.XYRange)
        xy = parameters.XYRange(i);
        measureLoop(xy,parameters)
    end % [end (xy) loop]
else
    for i = 1:length(parameters.XYRange)
        xy = parameters.XYRange(i);
        measureLoop(xy,parameters)
    end % [end (xy) loop]
end

% Before pool shutdown, pull number of workers (will be used to scale measurement combine)
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end

% Estimate total number of tracked cells
idx = 1;
while ~exist('CellData','var')
    i = parameters.XYRange(idx);
    parameters.XYDir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'xy',num2str(i),filesep]);
    if exist([parameters.XYDir,'CellData.mat'],'file')
        load([parameters.XYDir,'CellData.mat'])
    end
    idx = idx+1;
end
tot_cells = numel(parameters.XYRange)*size(CellData.FrameIn,1);


% Option 1: smaller experiment (50K cells) -> load ALL data and combine
if (poolsize==0) || (tot_cells<5e4)
    AllMeasurements = struct;
    AllMeasurements.parameters = orig_params;
    for i = parameters.XYRange
        parameters.XYDir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'xy',num2str(i),filesep]);
        if exist([parameters.XYDir,'CellMeasurements.mat'],'file')
            load([parameters.XYDir,'CellMeasurements.mat'])
            measureFields = fieldnames(CellMeasurements);
            for p = 1:length(measureFields)
                if (isfield(AllMeasurements,measureFields{p}))
                    AllMeasurements.(measureFields{p}) = cat(1, AllMeasurements.(measureFields{p}), CellMeasurements.(measureFields{p}));
                else
                    AllMeasurements.(measureFields{p}) = CellMeasurements.(measureFields{p});
                end
            end        
        end
    end % [end (xy) loop]
    save(namecheck([locations.data, filesep, parameters.SaveDirectory,filesep,'AllMeasurements.mat']),'AllMeasurements','-v7.3')
    
% Option 2: lg. experiment (>1e5 cells) -> load up data from (poolsize) cores @ a time, split fields into separate files
else
    boundary_pts = [1:poolsize:numel(parameters.XYRange), 1+numel(parameters.XYRange)];
    savedir = namecheck([locations.data, filesep, parameters.SaveDirectory,filesep,'AllMeasurements']);
    mkdir(savedir);
    for j = 1:(length(boundary_pts)-1)
        AllMeasurements = struct;
        % Part 1: combine data from (parpool size) sites
        for i = parameters.XYRange(boundary_pts(j):(boundary_pts(j+1)-1))
            parameters.XYDir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'xy',num2str(i),filesep]);
            if exist([parameters.XYDir,'CellMeasurements.mat'],'file')
                load([parameters.XYDir,'CellMeasurements.mat'])
                measureFields = fieldnames(CellMeasurements);
                for p = 1:length(measureFields)
                    if (isfield(AllMeasurements,measureFields{p}))
                        AllMeasurements.(measureFields{p}) = cat(1, AllMeasurements.(measureFields{p}), CellMeasurements.(measureFields{p}));
                    else
                        AllMeasurements.(measureFields{p}) = CellMeasurements.(measureFields{p});
                    end
                end        
            end
        end % [end (xy) loop]
            
            
        % Part II: load/create individual field files and write new data to them
        names = fieldnames(AllMeasurements);
        for i = 1:length(names)
            if j==1 
                if isnumeric(AllMeasurements.(names{i}))
                    eval([names{i}, '= single(AllMeasurements.(names{i}));'])
                else
                    eval([names{i}, '= AllMeasurements.(names{i});'])
                end
                save([savedir, filesep, names{i}, '.mat'],names{i},'-v7.3')
                clear(names{i});


            else % after 1st go-round, we need to append.
                load([savedir, filesep, names{i}, '.mat']);
                if isnumeric(AllMeasurements.(names{i}))
                    tmp =  single(AllMeasurements.(names{i}));
                else
                    tmp = AllMeasurements.(names{i});
                end
                eval([names{i},'= cat(1,',names{i},',tmp);']);  
                save([savedir, filesep, names{i}, '.mat'],names{i},'-v7.3')                         
                clear(names{i});
            end
        end
    end
    parameters = orig_params;
    save([savedir,filesep,'parameters.mat'],'parameters')
    
end
            
            
  
           
            
            


% Cycle through XY directories and combine their measurements



% Save AllMeasurements in condition directory. If we have trajectories for > 100K cells, save each field separately.
if size(AllMeasurements.CellData,1) < 1e5
    save(namecheck([locations.data, filesep, parameters.SaveDirectory,filesep,'AllMeasurements.mat']),'AllMeasurements','-v7.3')
else
    names = fieldnames(AllMeasurements);
    savedir = namecheck([locations.data, filesep, parameters.SaveDirectory,filesep,'AllMeasurements']);
    mkdir(savedir)
    for i =1:length(names)
        if isnumeric(AllMeasurements.(names{i}))
            eval([names{i}, '= single(AllMeasurements.(names{i}));'])
        else
            eval([names{i}, '= AllMeasurements.(names{i});'])
        end
        save([savedir, filesep, names{i}, '.mat'],names{i},'-v7.3')
        clear(names{i});
    end
end


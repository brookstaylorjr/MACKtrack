function [ScreenData, layout_dir, image_dir] = buildFromLayout(start_dir, nuclear_channel, measurement_channels, measurement_type, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [ScreenData, layout_dir, image_dir] = buildFromLayout(start_dir, nuclear_channel, measurement_channels, measurement_type, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% BUILDFROMLAYOUT crawls through all subdirectories under start_dir, looking for layout.xlsx files that correspond
% to plate screener images. Nuclear images will be segmented, and boundaries are used to measure intensity values
% (or other attribuites) of single cells in images specified in 'measurement_channels'
%
% A note about parameters - OP9_default works well for the 'standard' imaging conditions (20x, 2x2 binning, moderate 
% exposure), but you can modify/create new parameters by running the MACKtrack GUI
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Crawl down from starting directory - find all plate layouts (and corresponding image folders)
[layout_dir, image_dir, all_dir] = getImageDirectories(start_dir);

% For each layout, find all subsets (i.e. well(s) that correspond to conditions)
AllData = cell(size(layout_dir));
AllNames = cell(size(layout_dir));

for i = 1:length(layout_dir)
    [conditions, wells] = parseLayout(layout_dir{i});
    
    % DO MEASUREMENT HERE - assign to cell mat that matches layout_dir order (for assignment later)
    CondData = cell(size(conditions));
    CondName = cell(size(conditions));
    parfor k = 1:length(conditions)
        CondData{k} = segmentSubset(wells{k}, nuclear_channel, measurement_channels, ...
            measurement_type, image_dir{i}, parameters, [layout_dir{i},filesep,'output_img']);
        CondName{k} = conditions{k};
    end
    % Need to combine into a larger structure...
    AllData{i} = CondData;
    AllNames{i} = CondName;
end


% Write data into structure that matches folder hierarchy 
ScreenData= struct;
for i = 1:length(layout_dir)
    % Get structure field names (from file path)
    fieldnames = strsplit(layout_dir{i}((length(start_dir)+1):end),filesep);
    fieldnames(cellfun(@isempty,fieldnames)) = [];
    fieldnames = genvarname(fieldnames);
    final_name = 'SegmentData';
    % Make final structure
    if length(fieldnames)>0
        if ~isfield(ScreenData,fieldnames{1})
            ScreenData.(fieldnames{1}) = struct;
        end
        if length(fieldnames)>1
            for j = 2:length(fieldnames)
                structname = ['ScreenData.(''',strjoin(fieldnames(1:(j-1)),''').(''')];
                structname = [structname,''')'];
                if ~isfield(eval(structname),fieldnames{j})
                    eval([structname,'.(''',fieldnames{j},''') = struct;'])
                end
            end
        end
        final_name = ['ScreenData.(''',strjoin(fieldnames(1:(j)),''').('''), ''')'];
    
    end
    
    for j = 1:length(AllNames{i})
        eval([final_name,'.(''',AllNames{i}{j},''') = AllData{i}{j};'])
    
    end
    
end
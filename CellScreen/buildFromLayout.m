% Build measurment sets from layout files


% Specify staring directory - we'll find ALL layouts in any subdirectories under this one
start_dir = '/Volumes/labdata/prat';


% Specify channels and measurements
nuclear_channel = 'w1';
measurement_channels = {'w2','w5'};
measurement_type = {@mean,@median};


% Load parameters
load('OP9_prat.mat')

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
    parfor j = 1:length(conditions)
        CondData{j} = segmentSubset(wells{j}, nuclear_channel, measurement_channels, ...
            measurement_type, image_dir{i}, parameters, [layout_dir{i},filesep,'output_img']);
        CondName{j} = conditions{j};
    end
    % Need to combine into a larger structure...
    AllData{i} = CondData;
    AllNames{i} = CondName;
end

%%
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
        if ~isfield(fieldnames{1},ScreenData)
            ScreenData.(fieldnames{1}) = struct;
        end
        if length(fieldnames)>1
            for j = 2:length(fieldnames)
                structname = ['ScreenData.(''',strjoin(fieldnames(1:(j-1)),''').(''')];
                structname = [structname,''')'];
                if ~isfield(fieldnames{j},eval(structname))
                    eval([structname,'.(''',fieldnames{j},''') = struct;'])
                end
            end
        end
        final_name = ['ScreenData.(''',strjoin(fieldnames(1:(j-1)),''').('''), ''')'];
    
    end
    
    for j = 1:length(AllNames{i})
        eval([final_name,'.(''',AllNames{i}{j},''') = AllData{i}{j};'])
    
    end
    
end

save([start_dir,filesep,'ScreenData_',date,'.mat'],'ScreenData')


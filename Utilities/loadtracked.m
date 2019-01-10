function AllMeasurements = loadtracked(foldername, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% data = loadtracked(foldername, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LOADTRACKED pulls in AllMeasurements files (either the full file, or a few fields within
% a folder). By default, looks for a mounted external volume - otherwise, looks for server
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


external_name = '/Volumes/redstripe';
server_name = '/Volumes/labdata2/brooks';

    
stemlist = {['/Users/brooks/Desktop',filesep,foldername] %1st choice: check on desktop
    [external_name,filesep,'Tracked',filesep,foldername] % 2nd choice: look on external hard drive
        [server_name,filesep,'Tracked',filesep,foldername]
        ['Z:\brooks\Tracked',filesep,foldername]
        ['Y:\brooks\Tracked',filesep,foldername]
        }; % 3rd choice: look on server
   

%%
% Look on external HD 1st
for j = 1:length(stemlist)
    stem = stemlist{j};
    % 1st look for (smaller) live cell
    if exist([stem,filesep,'AllMeasurements.mat'],'file')
        load([stem,filesep,'AllMeasurements.mat'])
        return
    % 2nd look for fixed-cell
    elseif exist([stem,filesep,'AllData.mat'],'file')
        load([stem,filesep,'AllData.mat'])
        AllMeasurements = AllData;
        return
    % 3rd look for (larger) live cell    
    elseif exist([stem,filesep,'AllMeasurements'],'dir')
        if isempty(varargin)
            error('Please specify one or more fields (e.g. ''MeanNuc1'') to load')
        end
        AllMeasurements = struct;
        load([stem,filesep,'AllMeasurements',filesep,'CellData.mat']);
        AllMeasurements.CellData = double(CellData);
        for i = 1:length(varargin)
            tmp = load([stem,filesep,'AllMeasurements',filesep,varargin{i},'.mat']);
            if isnumeric(tmp.(varargin{i}))
                tmp.(varargin{i}) = double(tmp.(varargin{i}));
            end
            AllMeasurements = combinestructures(AllMeasurements,tmp); 
        end
        return
    end
end

error(['No ''',foldername, '''experiment folder found on server or external drive or desktop'])

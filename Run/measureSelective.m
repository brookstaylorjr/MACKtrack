measureSelective(parameters)


% Specify save directory and load spreadsheet URL
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end


% Grab sites that finished only!!
site_list = [];

for i = 1:length(parameters.XYRange)
    xyPos = parameters.XYRange(i);
    if exist(namecheck([locations.data,filesep, parameters.SaveDirectory,filesep,'xy',num2str(xyPos),filesep,'CellData.mat']),'file');
        site_list = cat(1,site_list,xyPos);
    end

end

disp(['measuring ',num2str(length(site_list),'sites (out of ',num2str(length(parameters.XYRange)) ' originally attempted)')

parameters.XYRange = site_list;

% Overwrite tracking params to only include sites that actually finished tracking
save(namecheck([locations.data,filesep, parameters.SaveDirectory,filesep,'TrackingParameters.mat']),'parameters')

% Measure selected sites
measureID(parameters)
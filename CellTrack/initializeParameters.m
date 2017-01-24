function [handlesOut] = initializeParameters(paramfile, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INITIALIZEPARAMETERS: update both handles structure GUI elements with loaded parameters
%
% paramfile  path to parameters file (string)
% handles   master structure with object handles and parameters
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

load(paramfile,'-mat')
handles.parameters = parameters;


% Load default parameters and fill in any fields that are missing in loaded file
load([handles.home_folder,'default_parameters.mat'],'-mat')
oldfields = fieldnames(parameters);
for i = 1:length(oldfields)
   if ~isfield(handles.parameters,oldfields{i})
        disp(['Creating field ',oldfields{i}])
        handles.parameters.(oldfields{i}) = parameters.(oldfields{i});
   end
end

% Loading/saving expressions and directories
set(handles.edit1A,'String',handles.parameters.ImagePath)
set(handles.edit2A,'String',handles.parameters.NucleusExpr)
set(handles.edit2B,'String',handles.parameters.CellExpr)
set(handles.edit2C,'String',handles.parameters.XYExpr)
set(handles.edit2D,'String',handles.parameters.TimeExpr)
set(handles.edit3A,'String',handles.parameters.SaveDirectory)

%% - - - - - - - - - Nuclear parameters - - - - - - - - - - - - 
nuc_edge = handles.parameters.NucleusEdgeThreshold;
set(handles.edit5A,'String',num2str(nuc_edge))
set(handles.slider5A,'Min',nuc_edge-1,'Max',nuc_edge+1, 'Value',nuc_edge) 
% Reset slider/axes ranges if necessary
if isfield(handles,'LineNuc1')
    x_lim = get(handles.axes5A,'XLim');
    if nuc_edge < x_lim(1)
        x_lim(1) = floor(nuc_edge);
    end
    if nuc_edge > x_lim(2)
        x_lim(2) = ceil(nuc_edge);
    end
    set(handles.slider5A,'Min',x_lim(1),'Max',x_lim(2)) 
    hold(handles.axes5A,'on')
    delete(handles.LineNuc1)
    handles.LineNuc1 = plot(handles.axes5A,ones(1,2)*nuc_edge,get(handles.axes5A,'YLim'),'Color',handles.blue);
    hold(handles.axes5A,'off')
    set(handles.axes5A,'XLim',x_lim);
end

set(handles.edit5B,'String',num2str(handles.parameters.WeakObjectCutoff))
set(handles.edit5C,'String',num2str(handles.parameters.MinNucleusRadius))
set(handles.edit5D,'String',num2str(handles.parameters.MaxNucleusRadius))
set(handles.popupmenu5A,'Value', find(strcmpi(get(handles.popupmenu5A,'String'),handles.parameters.ShapeDef),1,'first'));

set(handles.text5E,'String', [handles.parameters.ShapeDef,':']);
set(handles.edit5E,'String', num2str(handles.parameters.(handles.parameters.ShapeDef)(1)));
if length(handles.parameters.(handles.parameters.ShapeDef))<2
    val = handles.parameters.(handles.parameters.ShapeDef);
    handles.parameters.(handles.parameters.ShapeDef) = [val val];
end
set(handles.edit5F,'String', num2str(handles.parameters.(handles.parameters.ShapeDef)(2)));

set(handles.edit5G,'String',num2str(handles.parameters.NuclearSmooth))
set(handles.edit5H,'String',num2str(handles.parameters.NuclearInflection))


%% - - - - - - - - - Cell parameters - - - - - - - - - - - - 
set(handles.popupmenu0,'Value',find(strcmp(lower(handles.parameters.ImageType),{'phase','dic','fluorescence','none'})));
setVisibility(handles)

cell_edge1 = min(handles.parameters.CellSearchRange);
cell_edge2 = max(handles.parameters.CellSearchRange);
set(handles.edit6A,'String',num2str(cell_edge1))
set(handles.slider6A,'Min',cell_edge1-1,'Max',cell_edge2+1,'Value',cell_edge1)
set(handles.edit6B,'String',num2str(cell_edge2))
set(handles.slider6B,'Min',cell_edge1-1,'Max',cell_edge2+1,'Value',cell_edge2)

% Reset slider/axes ranges if necessary
if isfield(handles,'LineNoise1')
    x_lim = get(handles.axes6A,'XLim');
    if cell_edge1 < x_lim(1)
        x_lim(1) = floor(cell_edge1);
    end
    if cell_edge2 > x_lim(2)
        x_lim(2) = ceil(cell_edge2);
    end
    set(handles.slider6A,'Min',x_lim(1),'Max',x_lim(2))
    set(handles.slider6B,'Min',x_lim(1),'Max',x_lim(2)) 

    hold(handles.axes6A,'on')
    delete(handles.LineNoise1)
    delete(handles.LineNoise2)
    handles.LineNoise1 = plot(handles.axes6A,ones(1,2)*cell_edge1,get(handles.axes6A,'YLim'),'Color',handles.blue);
    handles.LineNoise2 = plot(handles.axes6A,ones(1,2)*cell_edge2,get(handles.axes6A,'YLim'),'Color',handles.blue);
    hold(handles.axes6A,'off')
    set(handles.axes6A,'XLim',x_lim);
end

set(handles.edit6C,'String',handles.parameters.MaxHoleSize)
set(handles.edit6D,'String',handles.parameters.HaloCutoff)
set(handles.edit6E,'String',handles.parameters.MinHoleSize)
set(handles.edit6F,'String',handles.parameters.MaxInflection)
set(handles.edit6G,'String',handles.parameters.WalkIn1)
set(handles.edit6H,'String',handles.parameters.WalkIn2)
set(handles.edit6I,'String',handles.parameters.LeakThrough)
gaussstring = '[';
for i = 1:length(handles.parameters.GaussSizes)
    gaussstring = [gaussstring,num2str(handles.parameters.GaussSizes(i)),','];
end
gaussstring(end) = ']';
set(handles.edit6J,'String',gaussstring)
set(handles.popupmenu6C,'Value',handles.parameters.Confluence+1)


%% - - - - - - - - - - Other parameters - - - - - - - - - - - - - -
set(handles.edit7A,'String',handles.parameters.StackSize)
set(handles.edit7B,'String',handles.parameters.DriftDistance)
set(handles.edit7D,'String',handles.parameters.MedianFilterSize)
set(handles.edit7E,'String',handles.parameters.NoiseSize)
set(handles.edit7G,'String',handles.parameters.MinCellWidth)
set(handles.edit7H,'String',handles.parameters.FramesPerHour)
jumpstring = '[';
for i = 1:length(handles.parameters.ImageJumps)
    jumpstring = [jumpstring,num2str(handles.parameters.ImageJumps(i)),','];
end
if strcmp(jumpstring(end),','); jumpstring = jumpstring(1:end-1); end
jumpstring = [jumpstring,']'];
set(handles.edit7J,'String',jumpstring)
% Get contents of 'CellMeasure' folder (the one stored with MACKtrack)
contents = dir([handles.home_folder,filesep, 'CellMeasure']);

% Find all 'modules' in directory, list
dropind = [];
for i = 1:length(contents)
    if isempty(strfind(lower(contents(i).name),'module')) % Capitalization-independent
        dropind = [dropind,i];
    else
        contents(i).name = contents(i).name(1:end-2); % (Drop the .m)
    end
end
contents(dropind) = [];
handles.parameters.ModuleNames = {contents.name};

% Create associated substructures/parameters (if required) for any discovered measurement modules
for i = 1:length(contents)
    if ~isfield(handles.parameters,contents(i).name)
       handles.parameters.(contents(i).name) = struct;
    end
end

check_fields = {'ImageExpr','ImageExpr2','ImageExpr3','Use'};
default_vals = {'--','--','--',get(handles.checkbox7A,'Min')};
for i = 1:length(contents)
    for j = 1:length(check_fields)
        if ~isfield(handles.parameters.(contents(i).name),check_fields{j})
            handles.parameters.(contents(i).name).(check_fields{j}) = default_vals{j};
        end
    end
end




% Cycle fields in parameters- if it refers to an old module, drop it.
testnames = fieldnames(handles.parameters);
for i = 1:length(testnames)
     if ~isempty(strfind(lower(testnames{i}),'module')) && ~strcmp(testnames{i},'ModuleNames')
         try 
             validatestring(testnames{i},handles.parameters.ModuleNames);
         catch me
            handles.parameters = rmfield(handles.parameters, testnames{i});
         end
     end
end

% Initialize/check measurement values in GUI 
set(handles.listbox7,'String',handles.parameters.ModuleNames)
set(handles.edit7K,'String',handles.parameters.(contents(1).name).ImageExpr)
set(handles.edit7L,'String',handles.parameters.(contents(1).name).ImageExpr2)
set(handles.edit7M,'String',handles.parameters.(contents(1).name).ImageExpr3)

set(handles.checkbox7A, 'Value', handles.parameters.(contents(1).name).Use);


% Set color of 'ok' based on string validity
try
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    filePath = [handles.locations.scope, handles.parameters.ImagePath,eval(handles.parameters.(contents(1).name).ImageExpr)];
    if exist(filePath,'file')
        set(handles.text7K_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7K_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7K_2,'ForegroundColor',handles.gray)
end

try
    filePath = [handles.locations.scope, handles.parameters.ImagePath,eval(handles.parameters.(contents(1).name).ImageExpr2)];
    if exist(filePath,'file')
        set(handles.text7L_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7L_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7L_2,'ForegroundColor',handles.gray)
end
try
    filePath = [handles.locations.scope, handles.parameters.ImagePath,eval(handles.parameters.(contents(1).name).ImageExpr3)];
    if exist(filePath,'file')
        set(handles.text7M_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7M_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7M_2,'ForegroundColor',handles.gray)
end


% Get flatfield objects and populate listbox accordingly
flatfields = {};
if isfield(handles.parameters,'Flatfield')
    for i = 1:length(handles.parameters.Flatfield)
        flatfields = cat(1,flatfields,{['Flatfield{',num2str(i),...
            '}(',num2str(size(handles.parameters.Flatfield{i},1)),' x ',...
            num2str(size(handles.parameters.Flatfield{i},1)),')']});
    end
end
set(handles.listbox7B,'String',flatfields)

% (Also update popupup menus)
if handles.parameters.CellFF > length(flatfields)
    handles.parameters.CellFF = 0;
    set(handles.popupmenu6A,'Value',0);
end
set(handles.popupmenu6B,'String',cat(1,{'None'},flatfields));
if handles.parameters.NucleusFF > length(flatfields)
    handles.parameters.NucleusFF = 0;
    set(handles.popupmenu6B,'Value',0);
end
set(handles.popupmenu6A,'String',cat(1,{'None'},flatfields));
set(handles.popupmenu6A,'Value',handles.parameters.NucleusFF+1)
set(handles.popupmenu6B,'Value',handles.parameters.CellFF+1)
guidata(handles.figure1,handles)





%% Resave parameters in place to reflect new updated values
load(paramfile,'-mat') % Reload original parameters for comparison
if ~isequal(handles.parameters,parameters)
    disp('Saving automatically-updated parameters')
    parameters = handles.parameters;
    save(paramfile,'parameters')
end


handlesOut = handles;

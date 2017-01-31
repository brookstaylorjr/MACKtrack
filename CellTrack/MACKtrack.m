function varargout = MACKtrack(varargin)
% MACKTRACK MATLAB code for MACKtrack.fig
%      MACKTRACK, by itself, creates a new MACKTRACK or raises the existing
%      singleton*.
%
%      H = MACKTRACK returns the handle to a new MACKTRACK or the handle to
%      the existing singleton*.
%
%      MACKTRACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MACKTRACK.M with the given input arguments.
%
%      MACKTRACK('Property','Value',...) creates a new MACKTRACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MACKtrack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MACKtrack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MACKtrack

% Last Modified by GUIDE v2.5 13-Jan-2017 16:52:28

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultuicontrolfontname','Sans Serif');
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MACKtrack_OpeningFcn, ...
                   'gui_OutputFcn',  @MACKtrack_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MACKtrack_OpeningFcn(hObject, eventdata, handles, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Executes just before MACKquery is made visible
%
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MACKquery (see VARARGIN)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Choose default command line output for MACKtrack
handles.output = hObject;

% Set colors
handles.blue = [0.051 0.329 0.671];
handles.orange = [207 79 51]/255;
handles.gray = [0.816 0.812 0.722];

% Specify MACKtrack home folder (one level above)
handles.home_folder = mfilename('fullpath');
slash_idx = strfind(handles.home_folder,filesep);
handles.home_folder = handles.home_folder(1:slash_idx(end-1));


% Set locations
load([handles.home_folder, 'locations.mat'],'-mat')

try
    locations.scope = namecheck(locations.scope,'');
    locations.data = namecheck(locations.data,'');
catch me
    h = setLocations;
    uiwait(h);
    load([handles.home_folder, 'locations.mat'],'-mat')
end

handles.locations = locations;
set(handles.text1A,'String',locations.scope);
set(handles.text3A,'String',locations.data);


% Set default parameters
handles = initializeParameters([handles.home_folder,'default_parameters.mat'], handles);


% Lock appropriate pushbuttons
handles.Locked = 1;
handles.Locked2 = 1;

% Load from parent directory and update handles structure
load_listbox(handles.parameters.ImagePath, handles)
handles = guidata(handles.figure1);
check_savedir(handles.parameters.SaveDirectory,handles);
handles = guidata(handles.figure1);

% OSX-specific formatting
% os = computer;
% if strcmp(os(1:4),'MACI')
%     a = fieldnames(handles);
%     for i = 1:length(a)
%         currentField = a{i};
%         if strcmp(currentField(1:4),'edit')
%             oldPos = get(handles.(currentField),'Position');
%             set(handles.(currentField),'Position', oldPos+[0 -1 0 2]);
%         end
%     end
% end
% ========================================================================================

function varargout = MACKtrack_OutputFcn(hObject, eventdata, handles) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set default command line output from handles structure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
varargout{1} = handles.output;
% ========================================================================================


%  - - - - - - - - - - - - - - - - UIPANEL 1 & 2 : FOLDER NAVIGATION/INPUT SELECTION  - - - - - - - -  - - - - - - - - - -

function edit1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT1A: Load listbox and update DataDirectory on selection
%
% hObject    handle to edit1A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

newfolder = get(hObject,'String');
if exist([handles.locations.scope,newfolder],'dir')
    % Remove leading slash, add trailing slash
    if  ~isempty(newfolder) && ~strcmp(newfolder(end),filesep)
        newfolder = [newfolder,filesep];
    end
    if ~isempty(newfolder) && strcmp(newfolder(1),filesep)
        newfolder = newfolder(2:end);
    end
    load_listbox(newfolder, handles)
else
    disp('- - - - - - - ')
    warning('Invalid folder entered')
    disp('- - - - - - - ')
end
% ========================================================================================

function listbox1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX1A: trigger directory update on double click/enter
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles  = guidata(handles.figure1);
% Determine user input
get(handles.figure1,'SelectionType');
% If user input is a double click, proceed
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox1A,'Value');
    file_list = get(handles.listbox1A,'String');
    % Get name of item selected in list box
    filename = file_list{index_selected};
    % If item is a directory, load list box with contents of new folder
    newfolder = [get(handles.edit1A,'String'),filename];
    % Remove leading slash, add trailing slash
    if  ~isempty(newfolder) && ~strcmp(newfolder(end),filesep)
        newfolder = [newfolder,filesep];
    end
    if ~isempty(newfolder) && strcmp(newfolder(1),filesep)
        newfolder = newfolder(2:end);
    end
    if exist([handles.locations.scope,newfolder],'dir')
        load_listbox(newfolder,handles)
    end
end
% ========================================================================================

function pushbutton1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1A: go up one directory level
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles  = guidata(handles.figure1);
newfolder = get(handles.edit1A,'String');

if ~isempty(newfolder) && strcmp(newfolder(end),filesep)
    newfolder =newfolder(1:end-1);
end
newfolder = newfolder(1:max(strfind(newfolder,filesep)));

% Remove leading slash, add trailing slash
if  ~isempty(newfolder) && ~strcmp(newfolder(end),filesep)
    newfolder = [newfolder,filesep];
end
if ~isempty(newfolder) && strcmp(newfolder(1),filesep)
    newfolder = newfolder(2:end);
end
load_listbox(newfolder,handles)

% ========================================================================================

function pushbutton1B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1B: browse for folder
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tmp = [handles.locations.scope,filesep,handles.parameters.ImagePath];
if ~exist(tmp,'dir')
    tmp = handles.locations.scope;
end
newfolder = uigetdir(tmp);
    
if (newfolder~=0)
    if length(newfolder)<length(handles.locations.scope) || ...
        ~strcmp(newfolder(1:length(handles.locations.scope)),handles.locations.scope)
        disp('- - - - - - - ')
        warning(['Images must be located at selected mount point (',handles.locations.scope,...
            ')- keeping previous directory.'])
        disp('- - - - - - - ')
    else
        newfolder = newfolder(length(handles.locations.scope)+1:end);
        % Remove leading slash, add trailing slash
        if  ~isempty(newfolder) && ~strcmp(newfolder(end),filesep)
            newfolder = [newfolder,filesep];
        end
        if ~isempty(newfolder) && strcmp(newfolder(1),filesep)
            newfolder = newfolder(2:end);
        end
        set(handles.edit1A,'String',newfolder)
        handles.parameters.ImagePath = newfolder;
        load_listbox(newfolder,handles)
    end
end
% ========================================================================================

function pushbutton1C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1C: specify new locations
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = setLocations;
uiwait(h);
load([handles.home_folder, 'locations.mat'],'-mat')
handles.locations = locations;
set(handles.text1A,'String',handles.locations.scope);
set(handles.text3A,'String',handles.locations.data);
drawnow;
load_listbox(handles.parameters.ImagePath,handles)
check_savedir(handles.parameters.SaveDirectory,handles)

% ========================================================================================


function load_listbox(dir_path, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LOAD_LISTBOX: Load listbox with contents of selected directory (subfunction)
%
% dir_path  name of directory to display
% handles   structure with handles and userdata
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ~exist([handles.locations.scope,dir_path],'dir')
    dir_path = '';
end

% Place remaining files in listbox and update edit1A (add filesep)
handles.file_names = quickdir([handles.locations.scope,dir_path]);

set(handles.listbox1A,'String',handles.file_names,...
 'Value',1)

set(handles.edit1A,'String',dir_path)
handles.parameters.ImagePath= dir_path;
check_expr(handles);
% ========================================================================================


function edit2A_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT2A-2D: set nuclear/cell img expressions, time and XY values, and update sample path
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
check_expr(handles);

function edit2B_Callback(hObject, eventdata, handles)
check_expr(handles);

function edit2C_Callback(hObject, eventdata, handles)
check_expr(handles);

function edit2D_Callback(hObject, eventdata, handles)
check_expr(handles);
% ========================================================================================


function check_expr(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECK_EXPR: Check/update paths to reflect image directory and name expression (subfunction)
%
% handles  structure with handles and user data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Don't let user proceed (or update parameter filename) unless file expression is valid
handles.Locked = 0;

% Ensure i/j vectors are valid
try
    % Define/check i and j indicies
    handles.parameters.XYExpr = get(handles.edit2C,'String');
    handles.parameters.TimeExpr = get(handles.edit2D,'String');    
   
    handles.parameters.XYRange = eval(handles.parameters.XYExpr);
    handles.parameters.TimeRange = eval(handles.parameters.TimeExpr);
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    
catch ME
    set(handles.text2G,'String','Error in i/j vector construction','ForegroundColor','r');
    set(handles.text2I,'String','Error in i/j vector construction','ForegroundColor','r');
    error('Invalid time or XY vector used')
    handles.Locked = 1;
end
    
   
% Check nuclear expression/ existence of sample nuclear image in parent directory
try 
    filestring = get(handles.edit2A,'String');
    sampleDirec = [handles.locations.scope, handles.parameters.ImagePath];
    sampleFile = eval(filestring);
    % Try to get an exact match (for tracking) or a partial one (for screens)
    pass = 1;
    partial_match = 0;
    if ~exist([handles.locations.scope,handles.parameters.ImagePath,sampleFile],'file')
        id = find(~cellfun(@isempty,strfind(handles.file_names,sampleFile)),1,'first');
        if ~isempty(id)
            id = find(~cellfun(@isempty,strfind(handles.file_names,sampleFile)));
            if length(id)>1
                id = id(cellfun(@isempty,strfind(handles.file_names(id),'thumb')));
                id = id(randperm(length(id),1));
                sampleFile = [handles.file_names{id}];
                partial_match = 1;
            end
            handles.parameters.NucleusMatch = filestring;
            filestring = ['''',sampleFile,''''];
        else
            pass = 0;
        end
    end
    if ~partial_match
        allowedLength = 75-length(sampleFile);
    else
        allowedLength = 64-length(sampleFile);
    end
    if length(sampleDirec)>allowedLength
        sampleDirec = [sampleDirec(1:floor(allowedLength/2 - 5)),'. . .',sampleDirec(end-ceil(allowedLength/2 - 5):end)];
    end
    if pass
        set(handles.text2G,'String',[sampleDirec,sampleFile],'ForegroundColor',handles.blue);
        handles.parameters.NucleusExpr = filestring;
        if partial_match
            set(handles.text2G,'String', [get(handles.text2G,'String'),' [EXAMP.]'])
        end 
    else
        set(handles.text2G,'String',['"',sampleFile,'" not found in current directory' ],'ForegroundColor','r');
        handles.Locked = 1;
    end

    

catch ME
     set(handles.text2G,'String','Invalid MATLAB string used for nuclear image expression','ForegroundColor','r');
     handles.Locked = 1;
end

% Check cell expression/ existence of sample cell file in parent directory (if applicable)
if ~strcmpi(handles.parameters.ImageType,'none')
    try 
        filestring = get(handles.edit2B,'String');
        sampleDirec = [handles.locations.scope, handles.parameters.ImagePath];
        sampleFile = eval(filestring);
        % Try to get an exact match (for tracking) or a partial one (for screens)
        pass = 1;
        if ~exist([handles.locations.scope,handles.parameters.ImagePath,sampleFile],'file')
            if ~isempty(id)
                id = find(~cellfun(@isempty,strfind(handles.file_names,sampleFile)));
                id = id(randperm(length(id),1));
                sampleFile = handles.file_names{id};
                handles.parameters.CellMatch = filestring;
                filestring = ['''',sampleFile,''''];
            else
                pass = 0;
            end
        end
    
        allowedLength = 75-length(sampleFile);
        if length(sampleDirec)>allowedLength
            sampleDirec = [sampleDirec(1:floor(allowedLength/2 - 5)),'. . .',sampleDirec(end-ceil(allowedLength/2 - 5):end)];
        end
        if pass
            set(handles.text2I,'String',[sampleDirec,sampleFile],'ForegroundColor',handles.blue);
            handles.parameters.CellExpr = filestring;        
        else
            set(handles.text2I,'String',['"',sampleFile,'" not found in current directory' ],'ForegroundColor','r');
            handles.Locked = 1;
        end
    catch ME
         set(handles.text2I,'String','Invalid MATLAB string used for cell image expression','ForegroundColor','r');
         handles.Locked = 1;
         ME
    end
else
    set(handles.text2I,'ForegroundColor',handles.gray);
    handles.Locked = 0;
end

% Change behavior based on lock
if handles.Locked
    set(handles.pushbutton2A,'ForegroundColor', handles.gray)
else
    set(handles.pushbutton2A,'ForegroundColor',handles.blue)
end
handles.Locked;
guidata(handles.figure1,handles)
% ========================================================================================


function partialmatch(string,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PARTIALMATCH if specified, then just matches first instance of string in a file list
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


function pushbutton2A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON2A: load images and create parameter visualizations
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ~handles.Locked
    set(handles.pushbutton2A,'ForegroundColor',handles.gray,'String','Loading...')
    drawnow;
    handles = loadImages(handles);
    handles.Locked2 = 0;
    set(handles.pushbutton2A,'ForegroundColor',handles.blue,'String','Load Images')
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
end
% ========================================================================================


%  - - - - - - - - - - - - - - - - UIPANEL 3: OUTPUT (SAVE) DIRECTORY - - - - - - - -  - - - - - - - - - -

function edit3A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT3A: select save directory for tracking output
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
check_savedir(get(handles.edit3A,'String'),handles);
% ========================================================================================

function pushbutton3A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON3A: Browse for folder
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tmp = [handles.locations.data,filesep,handles.parameters.SaveDirectory];
if ~exist(tmp,'dir')
    tmp = handles.locations.data;
end
newfolder = uigetdir(tmp);
if (newfolder~=0)
    if ~strcmp(newfolder(1:length(handles.locations.data)),handles.locations.data)
        warning(['Please choose a folder under the selected mount point (',handles.locations.data,')'])
    else
        check_savedir(newfolder(length(handles.locations.data)+1:end),handles);
    end
end

% ========================================================================================

function check_savedir(testfolder,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Helper function for UI panel 3: ensure that base location and subfolder location match.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
testfolder_orig = testfolder;
if ~isempty(testfolder) && strcmp(testfolder(end),filesep)
    testfolder = testfolder(1:end-1);
end
idx = strfind(testfolder,filesep);
if ~isempty(idx)
    idx= idx(end);
    testfolder = testfolder(1:idx);
else
    testfolder = '';
end
if ~exist([handles.locations.data, testfolder],'dir')
    if exist([handles.locations.data, handles.parameters.SaveDirectory],'dir')
        newfolder = handles.parameters.SaveDirectory;
        warning(['Data directory must be located under selected mount point (',handles.locations.data,...
            ')- resetting to previous SaveDirectory.'])
    else
        newfolder = '';
         warning(['Data directory must be located under selected mount point (',handles.locations.data,...
        ')- setting to base location.'])
        handles.parameters.SaveDirectory = '';
    end
else
    newfolder = testfolder_orig;

end

% Remove leading slash, add trailing slash
if  ~isempty(newfolder) && ~strcmp(newfolder(end),filesep)
    newfolder = [newfolder,filesep];
end
if ~isempty(newfolder) && strcmp(newfolder(1),filesep)
    newfolder = newfolder(2:end);
end
handles.parameters.SaveDirectory = newfolder;
set(handles.edit3A,'String',newfolder)
check_expr(handles)
guidata(handles.figure1,handles)

function pushbutton3B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON3B: specify new mount points for input/output data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = setLocations;
uiwait(h);
load([handles.home_folder, 'locations.mat'],'-mat')
handles.locations = locations;
set(handles.text1A,'String',locations.scope);
set(handles.text3A,'String',locations.data);
drawnow;
load_listbox(handles.parameters.ImagePath,handles)
check_savedir(handles.parameters.SaveDirectory,handles)
guidata(handles.figure1,handles)
% ========================================================================================



%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UIPANEL 4 : PARAMETERS AND TESTING/RUNNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pushbutton4A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON4A: Load an existing parameter set
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[FileName,PathName] = uigetfile('*','Select parameters file',[handles.home_folder,'Parameters']);
if (FileName~=0)
    handles = initializeParameters([PathName,FileName], handles);
end
% Load from parent directory and update handles structure
load_listbox(handles.parameters.ImagePath, handles)
% Resave parameters
parameters = handles.parameters;
save([PathName,FileName],'parameters');
% ========================================================================================

function pushbutton4B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON4B: save current parameter set
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
parameters = handles.parameters;
if ispc
    % Replace file separator in SaveDirectory and in ImagePath
    parameters.SaveDirectory(strfind(parameters.SaveDirectory,filesep)) = '/';
    parameters.ImagePath(strfind(parameters.ImagePath,filesep)) = '/';
end

FileName = '';
[FileName,PathName] = uiputfile([handles.home_folder,'Parameters',filesep,'*.mat'],...
    'Save current parameters (.mat)');
if (FileName~=0)
    save([PathName,FileName],'parameters')
end
% ========================================================================================

function pushbutton4C_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON4C: test first image in specified series
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ~handles.Locked2
    handles.parameters.debug = get(handles.checkbox4A,'Value');

    set(handles.pushbutton4C,'ForegroundColor',handles.gray,'String','Testing...')
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    drawnow;
    % Test images
    testImages(handles) % DIC or phase

    set(handles.pushbutton4C,'ForegroundColor',handles.blue,'String','Test')
    set(handles.pushbutton4D,'ForegroundColor',handles.blue, 'String','Run')
end
% ========================================================================================

function pushbutton4D_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON4D: run full analysis
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ~handles.Locked2
    try
        set(handles.pushbutton4C,'ForegroundColor',handles.gray)
        set(handles.pushbutton4D,'ForegroundColor',handles.gray,'String','Running...')
        drawnow;
        parameters = handles.parameters;    
        % 1) Track loop (put in parfor if 'parallel' option is selected'
        if get(handles.checkbox4B,'Value')
            parfor i = 1:length(parameters.XYRange)
                xyPos = parameters.XYRange(i);
                trackLoop(parameters,xyPos) % DIC or phase
            end        
            % 3) Measure loop (AllMeasurements.mat output)
            disp('Measuring...')
            MACKmeasure(parameters,1);

            set(handles.pushbutton4C,'ForegroundColor',handles.blue)
            set(handles.pushbutton4D,'ForegroundColor',handles.blue,'String','Run')

        else
            for i = 1:length(parameters.XYRange)
                xyPos = parameters.XYRange(i);
                trackLoop(parameters,xyPos) % DIC or phase
            end        
            % 3) Measure loop (AllMeasurements.mat output)
            disp('Measuring...')
            MACKmeasure(parameters,0);
        end

        set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue,'String','Run')
    catch ME
        set(handles.pushbutton4C,'ForegroundColor',handles.blue)
        set(handles.pushbutton4D,'ForegroundColor',handles.orange,'String','Run')
        rethrow(ME)
    end
end
% ========================================================================================

%  %%%%%%%%%%%%%%%%%% PARAMETER / IMAGE TYPE SWITCHING %%%%%%%%%%%%%%%%%%%%%%%%

function pushbuttonSwitch1_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTONSWITCH1: switch to nucleus parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set(handles.uipanel5,'Visible','on')
set(handles.uipanel6,'Visible','off')
set(handles.uipanel7,'Visible','off')

set(handles.pushbuttonSwitch1,'ForegroundColor','k')
set(handles.pushbuttonSwitch2,'ForegroundColor',handles.gray)
set(handles.pushbuttonSwitch3,'ForegroundColor',handles.gray)
% ========================================================================================

function pushbuttonSwitch2_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTONSWITCH2: switch to phase parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set(handles.uipanel5,'Visible','off')
set(handles.uipanel6,'Visible','on')
set(handles.uipanel7,'Visible','off')

set(handles.pushbuttonSwitch1,'ForegroundColor',handles.gray)
set(handles.pushbuttonSwitch2,'ForegroundColor','k')
set(handles.pushbuttonSwitch3,'ForegroundColor',handles.gray)
% ========================================================================================

function pushbuttonSwitch3_Callback(~, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTONSWITCH3: switch to segmentation/tracking parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set(handles.uipanel5,'Visible','off')
set(handles.uipanel6,'Visible','off')
set(handles.uipanel7,'Visible','on')
%set(handles.uipanel8,'Visible','off') (was image registration panel)

set(handles.pushbuttonSwitch1,'ForegroundColor',handles.gray)
set(handles.pushbuttonSwitch2,'ForegroundColor',handles.gray)
set(handles.pushbuttonSwitch3,'ForegroundColor','k')
% ========================================================================================


function popupmenu0_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POPUPMENU0: Switch between DIC and Phase-contrast algorithms and parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Get image type/update parameters
strlist = get(hObject,'String');
handles.parameters.ImageType = strlist{get(hObject,'Value')};
setVisibility(handles)
guidata(handles.figure1,handles)
% ========================================================================================





%  %%%%%%%%%%%%%%%%%%%     UIPANEL 5 : NUCLEUS PARAMETERS     %%%%%%%%%%%%%%%%%%%%%%


function slider5A_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SLIDER5A: set edge threshold for nucleus (NucleusEdgeThreshold)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
newVal = round(get(hObject,'Value'));
handles.parameters.NucleusEdgeThreshold = newVal;
set(handles.edit5A,'String',num2str(newVal))
set(handles.slider5A,'Value',newVal)
% Redo graph annotation
hold(handles.axes5A,'on')
 delete(handles.LineNuc1)
 handles.LineNuc1 = plot(handles.axes5A,ones(1,2)*newVal,get(handles.axes5A,'YLim'),'Color',handles.blue);
hold(handles.axes5A,'off')
guidata(handles.figure1,handles)
% ========================================================================================

function edit5A_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5A: set edge thresholding search range minimum (NucleusSearchRange- corresponds with slider5A)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
try
    newVal = round(eval(get(hObject,'String')));
    if newVal > get(handles.slider5A,'Max')
        newVal = get(handles.slider5A,'Max');
    end
    if newVal < get(handles.slider5A,'Min')
        newVal = get(handles.slider5A,'Min');
    end

    handles.parameters.NucleusEdgeThreshold = newVal;
    set(handles.edit5A,'String',num2str(newVal))
    set(handles.slider5A,'Value',newVal)
    % Redo graph annotation
    hold(handles.axes5A,'on')
    delete(handles.LineNuc1)
    handles.LineNuc1 = plot(handles.axes5A,ones(1,2)*newVal,get(handles.axes5A,'YLim'),'Color',handles.blue);
    hold(handles.axes5A,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)

catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end



function edit5C_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5C: set nuclear size limit (MinNucleusRadius- displays on axes5B)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MinNucleusRadius = newVal;
    nucCircle1 = getnhood(diskstrel(handles.parameters.MinNucleusRadius));
    nucCircle2 = getnhood(diskstrel(handles.parameters.MaxNucleusRadius));
    diskNuc = padarray(nucCircle1,floor((handles.parameters.ImageSize-size(nucCircle1))/2),0,'pre' );
    diskNuc_temp = padarray(nucCircle2,floor((handles.parameters.ImageSize-size(nucCircle2))/2),0,'pre' );
    diskNuc = padarray(diskNuc,ceil((handles.parameters.ImageSize-size(nucCircle1))/2),0,'post' ) + ...
    padarray(diskNuc_temp,ceil((handles.parameters.ImageSize-size(nucCircle2))/2),0,'post' );
    set(handles.figure1,'CurrentAxes',handles.axes5B)
    hold(handles.axes5B,'on')
    delete(handles.NucleusOverlay)
    handles.NucleusOverlay = imshow(label2rgb(diskNuc,'spring')); 
    set(handles.NucleusOverlay,'AlphaData',double(diskNuc>0)*0.6)
    hold(handles.axes5B,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit5B_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5B: set weak object morphological cutoff (WeakObjectCutoff)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

try
    newVal = eval(get(hObject,'String'));
    set(hObject,'String',num2str(newVal))
    handles.parameters.WeakObjectCutoff = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================



function edit5D_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5D: set nuclear size limit (MaxNucleusRadius- displays on axes5B)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MaxNucleusRadius = newVal;
    nucCircle1 = getnhood(diskstrel(handles.parameters.MinNucleusRadius));
    nucCircle2 = getnhood(diskstrel(handles.parameters.MaxNucleusRadius));
    diskNuc = padarray(nucCircle1,floor((handles.parameters.ImageSize-size(nucCircle1))/2),0,'pre' );
    diskNuc_temp = padarray(nucCircle2,floor((handles.parameters.ImageSize-size(nucCircle2))/2),0,'pre' );
    diskNuc = padarray(diskNuc,ceil((handles.parameters.ImageSize-size(nucCircle1))/2),0,'post' ) + ...
    padarray(diskNuc_temp,ceil((handles.parameters.ImageSize-size(nucCircle2))/2),0,'post' );
    set(handles.figure1,'CurrentAxes',handles.axes5B)
    hold(handles.axes5B,'on')
    delete(handles.NucleusOverlay)
    handles.NucleusOverlay = imshow(label2rgb(diskNuc,'spring')); 
    set(handles.NucleusOverlay,'AlphaData',double(diskNuc>0)*0.6)
    hold(handles.axes5B,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function popupmenu5A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POPUP5A: choose parameter to define nuclear shape: compactness or solidity
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
idx = get(hObject,'Value');
fullstr = get(hObject,'String');
handles.parameters.ShapeDef = fullstr{idx};
set(handles.text5E,'String', [handles.parameters.ShapeDef,':']);
set(handles.edit5E,'String', num2str(handles.parameters.(handles.parameters.ShapeDef)(1)));
try
    set(handles.edit5F,'String', num2str(handles.parameters.(handles.parameters.ShapeDef)(2)));
catch me
    set(handles.edit5F,'String', num2str(handles.parameters.(handles.parameters.ShapeDef)(1)));
    handles.parameters.(handles.parameters.ShapeDef)(2) = handles.parameters.(handles.parameters.ShapeDef)(1);
end
guidata(handles.figure1,handles)



function edit5E_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5E: set strict compactness/solidity bound (for round nuclei)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    if strcmp(handles.parameters.ShapeDef,'Compactness')
        test1 = @(v1,v2) v1>v2;
    else
        test1 = @(v1,v2) v1<v2;
    end
    if test1(newVal,handles.parameters.(handles.parameters.ShapeDef)(2))
        warning('Resetting value - lenient compactness should be higher than strict, lenient soliditiy should be lower')
        newVal = handles.parameters.(handles.parameters.ShapeDef)(2);
    end
    set(hObject,'String',num2str(newVal))
    handles.parameters.(handles.parameters.ShapeDef)(1) = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit5F_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5F: set lenient compactness/solidity bound (for oblong nuclei)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    % Parameter checking
    if strcmp(handles.parameters.ShapeDef,'Compactness')
        test1 = @(v1,v2) v1<v2;
    else
        test1 = @(v1,v2) v1>v2;
    end
    if test1(newVal,handles.parameters.(handles.parameters.ShapeDef)(1))
        warning('Resetting value - lenient compactness should be higher than strict, lenient soliditiy should be lower')
        newVal = handles.parameters.(handles.parameters.ShapeDef)(1);
    end
    set(hObject,'String',num2str(newVal))
    handles.parameters.(handles.parameters.ShapeDef)(2) = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit5G_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5G: set nuclear smoothing size (nan is ok)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    % Parameter checking
    if ~isnumeric(newVal)
        warning('Resetting - value must be numeric')
        newVal = handles.parameters.NuclearSmooth;
    end
    set(hObject,'String',num2str(newVal))
    handles.parameters.NuclearSmooth = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================



function edit5H_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT5H: set minimum inflection angle for splitting nuclei
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    % Parameter checking
    if ~isnumeric(newVal)
        warning('Resetting - value must be numeric')
        newVal = handles.parameters.NuclearInflection;
    end
    set(hObject,'String',num2str(newVal))
    handles.parameters.NuclearInflection = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================



% %%%%%%%%%%%%%%%%%%%%%% UIPANEL 6 : Phase/DIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slider6A_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SLIDER6A: set noise thresholding search range minimum (CellSearchRange)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
newVal = round(get(hObject,'Value'));
if newVal > get(handles.slider6B,'Value')
    newVal = get(handles.slider6B,'Value');
end
handles.parameters.CellSearchRange(1) = newVal;
set(handles.edit6A,'String',num2str(newVal))
set(handles.slider6A,'Value',newVal)
% Redo graph annotation
hold(handles.axes6A,'on')
 delete(handles.LineNoise1)
 handles.LineNoise1 = plot(handles.axes6A,ones(1,2)*newVal,get(handles.axes6A,'YLim'),'Color',handles.blue);
hold(handles.axes6A,'off')
guidata(handles.figure1,handles)
% ========================================================================================

function edit6A_Callback(hObject, ~, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6A: set noise thresholding search range minimum (CellSearchRange- corresponds with slider5A)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
try
    newVal = round(eval(get(hObject,'String')));
    if newVal > get(handles.slider6A,'Max')
        newVal = get(handles.slider6A,'Max');
    end
    if newVal < get(handles.slider6A,'Min')
        newVal = get(handles.slider6A,'Min');
    end
    if newVal > get(handles.slider6B,'Value')
        newVal = get(handles.slider6B,'Value');
    end
    handles.parameters.CellSearchRange(1) = newVal;
    set(handles.edit6A,'String',num2str(newVal))
    set(handles.slider6A,'Value',newVal)
    % Redo graph annotation
    hold(handles.axes6A,'on')
     delete(handles.LineNoise1)
     handles.LineNoise1 = plot(handles.axes6A,ones(1,2)*newVal,get(handles.axes6A,'YLim'),'Color',handles.blue);
    hold(handles.axes6A,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)

catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function slider6B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SLIDER6B: set noise thresholding search range maximum (CellSearchRange)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
newVal = round(get(hObject,'Value'));
if newVal < get(handles.slider6A,'Value')
    newVal = get(handles.slider6A,'Value');
end
handles.parameters.CellSearchRange(2) = newVal;
set(handles.edit6B,'String',num2str(newVal))
set(handles.slider6B,'Value',newVal)
% Redo graph annotation
hold(handles.axes6A,'on')
 delete(handles.LineNoise2)
 handles.LineNoise2 = plot(handles.axes6A,ones(1,2)*newVal,get(handles.axes6A,'YLim'),'Color',handles.blue);
hold(handles.axes6A,'off')
guidata(handles.figure1,handles)
% ========================================================================================

function edit6B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6B: set noise thresholding search range minimum (CellSearchRange- corresponds with slider5A)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get value/ update GUI elements
try
    newVal = round(eval(get(hObject,'String')));
    if newVal > get(handles.slider6A,'Max')
        newVal = get(handles.slider6A,'Max');
    end
    if newVal < get(handles.slider6A,'Min')
        newVal = get(handles.slider6A,'Min');
    end
    if newVal < get(handles.slider6A,'Value')
        newVal = get(handles.slider6A,'Value');
    end
    handles.parameters.CellSearchRange(2) = newVal;
    set(handles.edit6B,'String',num2str(newVal))
    set(handles.slider6B,'Value',newVal)
    % Redo graph annotation
    set(handles.figure1,'CurrentAxes',handles.axes6A)
    hold(handles.axes6A,'on')
     delete(handles.LineNoise2)
     handles.LineNoise2 = plot(handles.axes6A,ones(1,2)*newVal,get(handles.axes6A,'YLim'),'Color',handles.blue);
    hold(handles.axes6A,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)

catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit6C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6C: set max filled hole size (MaxHoleSize - displays on axes6B)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MaxHoleSize = newVal;
    % Plot annotation disk: CellOverlay1/CellOverlay2
    parameters = handles.parameters;
    radius1 = round(sqrt(parameters.MaxHoleSize/pi));
    cellCircle1 = getnhood(diskstrel(radius1));
    radius2 = round(sqrt(parameters.MinHoleSize/pi));
    cellCircle2 = getnhood(diskstrel(radius2));
    diskCell = padarray(cellCircle1,floor((parameters.ImageSize-size(cellCircle1))/2),0,'pre' );
    diskCell_temp = padarray(cellCircle2,floor((parameters.ImageSize-size(cellCircle2))/2),0,'pre' );
    diskCell = padarray(diskCell,ceil((parameters.ImageSize-size(cellCircle1))/2),0,'post' ) + ...
        padarray(diskCell_temp,ceil((parameters.ImageSize-size(cellCircle2))/2),0,'post' );
    set(handles.figure1,'CurrentAxes',handles.axes6B)
    hold(handles.axes6B,'on')
    delete(handles.CellOverlay)
    handles.CellOverlay = imshow(label2rgb(diskCell,'winter'),'Parent',handles.axes6B); 
    set(handles.CellOverlay,'AlphaData',double(diskCell>0)*0.3)
    hold(handles.axes6B,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit6D_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6D: set cutoff proportion of thresholded "bright" pixels (HaloCutoff)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.HaloCutoff = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit6E_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6E: set minimum filled hole size (MinHoleSize - displays on axes6B)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MinHoleSize = newVal;
    % Plot annotation disk: CellOverlay1/CellOverlay2
    parameters = handles.parameters;
    radius1 = round(sqrt(parameters.MaxHoleSize/pi));
    cellCircle1 = getnhood(diskstrel(radius1));
    radius2 = round(sqrt(parameters.MinHoleSize/pi));
    cellCircle2 = getnhood(diskstrel(radius2));
    diskCell = padarray(cellCircle1,floor((parameters.ImageSize-size(cellCircle1))/2),0,'pre' );
    diskCell_temp = padarray(cellCircle2,floor((parameters.ImageSize-size(cellCircle2))/2),0,'pre' );
    diskCell = padarray(diskCell,ceil((parameters.ImageSize-size(cellCircle1))/2),0,'post' ) + ...
        padarray(diskCell_temp,ceil((parameters.ImageSize-size(cellCircle2))/2),0,'post' );
    set(handles.figure1,'CurrentAxes',handles.axes6B)
    hold(handles.axes6B,'on')
    delete(handles.CellOverlay)
    handles.CellOverlay = imshow(label2rgb(diskCell,'winter'),'Parent',handles.axes6B); 
    set(handles.CellOverlay,'AlphaData',double(diskCell>0)*0.3)
    hold(handles.axes6B,'off')
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit6F_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6F: set maximum distance that inflection-based segmentation is applied (MaxInflection)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MaxInflection = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================



function edit6G_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6G: set number of steps for 1st stage of walk-in (WalkIn1)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.WalkIn1 = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit6H_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6H: set number of steps for 2nd stage of walk-in (WalkIn2)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.WalkIn2 = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit6I_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6I: set minimum size for walk-in edges (LeakThrough- contiguous pixels)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.LeakThrough = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit6J_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6J: set vector of Gaussian radii for scale space edge-finding 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    gaussstring = '[';
    for i = 1:length(newVal)
        gaussstring = [gaussstring,num2str(newVal(i)),','];
    end
    gaussstring(end) = ']';
    set(hObject,'String',gaussstring)
    handles.parameters.GaussSizes = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end

% ========================================================================================

function popupmenu6A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POPUPMENU6A: specify flatfield correction image for nuclear channel (NucleusFF)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
handles.parameters.NucleusFF = get(hObject,'Value') - 1;
guidata(handles.figure1,handles)
% ========================================================================================


function popupmenu6B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POPUPMENU6B: specify flatfield correction image for cell channel (CellFF)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
handles.parameters.CellFF = get(hObject,'Value') - 1;
guidata(handles.figure1,handles)
% ========================================================================================

function popupmenu6C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POPUPMENU6C: specify image type (used to choose a thresholding method) (Confluence)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
handles.parameters.Confluence = get(hObject,'Value') - 1;
guidata(handles.figure1,handles)
% ========================================================================================


%  - - - - - - - - - - - - - - -  - - - - - - UIPANEL 7 : OTHER PARAMETERS - - - - - - - - - - - - - - - - - - - - - - -

function edit7A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7A: set stack size for tracking process (StackSize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.StackSize = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit7B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7B: set maximum frame-to-frame drift distance (DriftDistance)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.DriftDistance = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit7D_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7D: square median filter size (MedianFilterSize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MedianFilterSize = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit7E_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7E: image noise size (NoiseSize- contiguous pixels)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.NoiseSize = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================

function edit7G_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7G: set minumum cell width (MinCellWidth- round strel diameter)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = round(eval(get(hObject,'String')));
    set(hObject,'String',num2str(newVal))
    handles.parameters.MinCellWidth = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function edit7H_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7H: set frames imaged per hour - default is 12 (FramesPerHour)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    set(hObject,'String',num2str(newVal))
    handles.parameters.FramesPerHour = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================



function edit7J_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7J: vector of timepoints where image "jumped" (activates cross-correlation normalization 
% for these images) (ImageJumps)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newVal = eval(get(hObject,'String'));
    jumpstring = '[';
    for i = 1:length(newVal)
        jumpstring = [jumpstring,num2str(newVal(i)),','];
    end
    if strcmp(jumpstring(end),','); jumpstring = jumpstring(1:end-1); end
    jumpstring = [jumpstring,']'];    
    set(hObject,'String',jumpstring)
    handles.parameters.ImageJumps = newVal;
    handles.Locked2 = 0;
    set(handles.pushbutton4C,'ForegroundColor',handles.blue)
    set(handles.pushbutton4D,'ForegroundColor',handles.blue)
    guidata(handles.figure1,handles)
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.pushbutton4C,'ForegroundColor',handles.gray)
    set(handles.pushbutton4D,'ForegroundColor',handles.gray)
    rethrow(ME)
end
% ========================================================================================


function listbox7_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX7: Load module's associated parameters and update GUI elements
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get selected module
index_selected = get(handles.listbox7,'Value');
module_list = get(handles.listbox7,'String');
module = module_list{index_selected};

set(handles.edit7K,'String',handles.parameters.(module).ImageExpr);


try
    set(handles.edit7L,'String',handles.parameters.(module).ImageExpr2);
catch ME
    handles.parameters.(module).ImageExpr2 = '--';
    guidata(handles.figure1,handles)
    set(handles.edit7L,'String',handles.parameters.(module).ImageExpr2);
end

set(handles.checkbox7A, 'Value', handles.parameters.(module).Use);

% Set color of 'ok' based on string validity
try
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    filePath = [handles.locations.scope,handles.parameters.ImagePath,eval(handles.parameters.(module).ImageExpr)];
    if exist(filePath,'file')
        set(handles.text7K_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7K_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7K_2,'ForegroundColor',handles.gray)
end
try
    filePath = [handles.locations.scope,handles.parameters.ImagePath,eval(handles.parameters.(module).ImageExpr2)];
    if exist(filePath,'file')
        set(handles.text7L_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7L_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7L_2,'ForegroundColor',handles.gray)
end
% ========================================================================================

function checkbox7A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX7A: check/set measurement use for selected module
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
% Get selected module
index_selected = get(handles.listbox7,'Value');
module_list = get(handles.listbox7,'String');
module = module_list{index_selected};
 
handles.parameters.(module).Use = get(hObject,'Value');
guidata(handles.figure1,handles)


function edit7K_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7K: check/set measurement image expression for selected module
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get selected module
index_selected = get(handles.listbox7,'Value');
module_list = get(handles.listbox7,'String');
module = module_list{index_selected};
handles.parameters.(module).ImageExpr = get(hObject,'String');   
guidata(handles.figure1,handles)

try
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    filePath = [handles.locations.scope,handles.parameters.ImagePath,eval(get(hObject,'String'))];
    if exist(filePath,'file')
        set(handles.text7K_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7K_2,'ForegroundColor',handles.gray)
    end
catch ME
    set(handles.text7K_2,'ForegroundColor',handles.gray)
    disp(ME.message)
    
end
% ========================================================================================



function edit7L_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7L: check/set secondary measurement image expression for selected module
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get selected module
index_selected = get(handles.listbox7,'Value');
module_list = get(handles.listbox7,'String');
module = module_list{index_selected};
handles.parameters.(module).ImageExpr2 = get(hObject,'String');   
guidata(handles.figure1,handles)
try
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    filePath = [handles.locations.scope,handles.parameters.ImagePath,eval(get(hObject,'String'))];
    if exist(filePath,'file')
        set(handles.text7L_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7L_2,'ForegroundColor',handles.gray)
    end
    
catch ME
    set(handles.text7L_2,'ForegroundColor',handles.gray)
    disp(ME.message)
end
% ========================================================================================

function edit7M_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7M: check/set tertiary measurement image expression for selected module
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get selected module
index_selected = get(handles.listbox7,'Value');
module_list = get(handles.listbox7,'String');
module = module_list{index_selected};
handles.parameters.(module).ImageExpr3 = get(hObject,'String');   
guidata(handles.figure1,handles)
try
    i = min(handles.parameters.XYRange);
    j = min(handles.parameters.TimeRange);
    filePath = [handles.locations.scope,handles.parameters.ImagePath,eval(get(hObject,'String'))];
    if exist(filePath,'file')
        set(handles.text7M_2,'ForegroundColor',handles.blue)    
    else
        set(handles.text7M_2,'ForegroundColor',handles.gray)
    end
    
catch ME
    set(handles.text7M_2,'ForegroundColor',handles.gray)
    disp(ME.message)
end


function listbox7B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX7A: show all added flatfield images - show image/image location if item is double-clicked
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
get(handles.figure1,'SelectionType');
% If user input is a double click, proceed
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox7B,'Value');
    figure,imagesc(handles.parameters.Flatfield{index_selected})
    set(gca,'XTick',[],'YTick',[]), colorbar
    axis equal
    title1 = handles.parameters.FlatfieldNames{index_selected};
    title2 = {};
    splits = [1:96:length(title1),length(title1)+1];
    for i = 1:(length(splits)-1)
    title2 = cat(1,title2,{title1(splits(i):(splits(i+1)-1))});
    end
    title(title2,'Interpreter','none',...
        'FontSize',8,'FontWeight','normal')
end


function pushbutton7A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON7A: add a new flatfield image (stored by name and image)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
[file1, dir1] = uigetfile({'*.tif;*.tiff;*.png;*.jpg','Image Files';'*.*', 'All Files'},'Load a flatfield image');
if file1
    % Add file/filename to parameters
    new_img = double(imread([dir1,filesep,file1]));
    if ~isfield(handles.parameters,'Flatfield')
        handles.parameters.Flatfield{1} = new_img;
        handles.parameters.FlatfieldNames{1} = [dir1,filesep,file1];
        set(handles.listbox7B,'Value',1)
    else
        handles.parameters.Flatfield = cat(2,handles.parameters.Flatfield, {new_img});
        handles.parameters.FlatfieldNames = cat(2,handles.parameters.FlatfieldNames,{[dir1,filesep,file1]});
    end
    updateflatfields(handles);
    guidata(handles.figure1,handles);
end

function pushbutton7B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON7B: remove selected flatfield image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
index_selected = get(handles.listbox7B,'Value');
handles.parameters.Flatfield(index_selected) = [];
handles.parameters.FlatfieldNames(index_selected) = [];
disp('Deleted uploaded flatfield image')
updateflatfields(handles);
guidata(handles.figure1,handles);



function pushbutton7C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON7C: swap selected flatfield image for a new image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
index_selected = get(handles.listbox7B,'Value');
[file1, dir1] = uigetfile({'*.tif;*.tiff;*.png;*.jpg','Image Files';'*.*', 'All Files'},'Load a flatfield image',...
    [handles.locations.scope,handles.parameters.ImagePath]);
if file1
    % Add file/filename to parameters
    new_img = double(imread([dir1,filesep,file1]));
    handles.parameters.Flatfield{index_selected} = new_img;
    handles.parameters.FlatfieldNames{index_selected} = [dir1,filesep,file1];
    % Update listbox and save parameters
    updateflatfields(handles);
end


function updateflatfields(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Update list of flatfields based on parameters (similar to lines in initializeParameters.m)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
flatfields = {};
if isfield(handles.parameters,'Flatfield')
    for i = 1:length(handles.parameters.Flatfield)
        flatfields = cat(1,flatfields,{['Flatfield{',num2str(i),...
            '}(',num2str(size(handles.parameters.Flatfield{i},1)),' x ',...
            num2str(size(handles.parameters.Flatfield{i},1)),')']});
    end
end
set(handles.listbox7B,'String',flatfields)
% Also update popupup menus accordingly
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
guidata(handles.figure1,handles)


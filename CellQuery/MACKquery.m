function varargout = MACKquery(varargin)
% MACKquery M-file for MACKquery.fig
%      MACKQUERY, by itself, creates a new MACKQUERY or raises the existing
%      singleton*.
%
%      H = MACKQUERY returns the handle to a new MACKQUERY or the handle to
%      the existing singleton*.
%
%      MACKQUERY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MACKQUERY.M with the given input arguments.
%
%      MACKQUERY('Property','Value',...) creates a new MACKQUERY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MACKquery_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MACKquery_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MACKquery

% Last Modified by GUIDE v2.5 10-Apr-2015 13:43:06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin initialization code - DO NOT EDIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MACKquery_OpeningFcn, ...
                   'gui_OutputFcn',  @MACKquery_OutputFcn, ...
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

% ______________________________________________ OPENING/OUTPUT FUNCTIONS ______________________________________________

function MACKquery_OpeningFcn(hObject, eventdata, handles, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Executes just before MACKquery is made visible
%
% This function has no output args, see ExportFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MACKquery (see VARARGIN)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Choose default command line output for MACKquery
handles.output = hObject;

% Load up default locations
handles.home_folder = mfilename('fullpath');
slash_idx = strfind(handles.home_folder,filesep);
handles.home_folder = handles.home_folder(1:slash_idx(end-1));

% Set locations
load([handles.home_folder, 'locations.mat'],'-mat')
handles.locations = locations;

% Set default parameters and load from starting directory
handles = setOptions(handles);
guidata(handles.figure1,handles)
load_listbox(handles.Options.DataDirectory, handles)

% Mac-specific formatting:
os = computer;
if strcmp(os(1:4),'MACI')
a = fieldnames(handles);
for i = 1:length(a)
    currentField = a{i};
    if strcmp(currentField(1:4),'edit')
        oldPos = get(handles.(currentField),'Position');
        set(handles.(currentField),'Position', oldPos+[0 -3 0 3]);
    end
end
end


% ========================================================================================

function varargout = MACKquery_OutputFcn(hObject, eventdata, handles) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set default command line output from handles structure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
varargout{1} = handles.output;
% ========================================================================================




%  - - - - - - - - - - - - - - - - - - - - UIPANEL 1 & 2: FOLDER NAVIGATION  - - - - - - - - - - - - - - - - - - - -

function edit1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT1A: Load listbox and update DataDirectory on selection
%
% hObject    handle to edit1A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.Options.DataDirectory = get(hObject,'String');
if exist(handles.Options.DataDirectory,'dir')
    load_listbox(handles.Options.DataDirectory, handles)
end
% ========================================================================================

function listbox1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX1A: trigger directory update on double click/enter
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% First determine user input
get(handles.figure1,'SelectionType');
% If user input is a double click
if strcmp(get(handles.figure1,'SelectionType'),'open')
    index_selected = get(handles.listbox1A,'Value');
    file_list = get(handles.listbox1A,'String');
    % Get name of item selected in list box
    filename = file_list{index_selected};
    % If item is a directory, load list box with contents of new folder
    if  handles.is_dir(handles.sorted_index(index_selected))
        newFolder = get(handles.edit1A,'String');
        if ~strcmp(newFolder(end),filesep)
            newFolder = [newFolder,filesep];
        end
        newFolder = [newFolder,filename];
        load_listbox(newFolder,handles)
    % Alternately, if AllMeasurements.mat is double clicked, load measurements
    elseif ~isempty(strfind(filename,'.mat'))
        load_measurements(filename, handles);
    end
end
% ========================================================================================

function pushbutton1A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1A: go up one directory level
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

newFolder = get(handles.edit1A,'String'); 
if strcmp(newFolder(end),filesep)
            newFolder =newFolder(1:end-1);
end
newFolder = newFolder(1:max(strfind(newFolder,filesep)));
load_listbox(newFolder,handles)

% ========================================================================================

function pushbutton1B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1B: save parent directory and measurements, update "Contents" text. Ensure 
% selected file is "AllMeasurements.mat"- if not, display warning.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Get name of item selected in list box
index_selected = get(handles.listbox1A,'Value');
file_list = get(handles.listbox1A,'String');
filename = file_list{index_selected};
if ~isempty(strfind(filename,'.mat'))
    load_measurements(filename,handles);
else
    disp('Maybe try to pick a .mat file?')
end

% ========================================================================================

function pushbutton1C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1C: browse for folder
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
PathName = uigetdir(get(handles.edit1A,'String'));
if ~strcmp(PathName(end),filesep)
	PathName = [PathName,filesep];
end
if (PathName~=0)
    set(handles.edit1A,'String',PathName)
    handles.Options.DataDirectory = PathName;
    load_listbox(PathName,handles)
end
% ========================================================================================

function load_listbox(dir_path, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Load listbox with contents of selected directory (subfunction)
%
% dir_path  name of directory to display
% handles   structure with handles and userdata
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

set(handles.listbox1A,'String','Loading...')
dir_struct = dir(dir_path);
% Clean directory of any invisible files
drop_ind = 1;
for i=1:length(dir_struct)
    if strcmp(dir_struct(drop_ind).name(1),'.')
        dir_struct(drop_ind) = [];
    else
        drop_ind = drop_ind+1;
    end
end
% Place remaining files in listbox and update edit1A
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.figure1,handles)
set(handles.listbox1A,'String',handles.file_names,...
 'Value',1)
set(handles.edit1A,'String',dir_path)
handles.Options.DataDirectory = get(handles.edit1A,'String');
guidata(handles.figure1,handles)
% ========================================================================================

function load_measurements(filename, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LOADMEASUREMENTS checks .mat file is valid, and loads measurments into handles, pulling
% off parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
load([handles.Options.DataDirectory,filesep,filename]);
if exist('AllMeasurements','var')
    handles.file_loc = [handles.Options.DataDirectory,filesep,filename];
    % Pull off tracking parameters and save
    if isfield(AllMeasurements,'parameters')
        handles.trackparameters = AllMeasurements.parameters;
        AllMeasurements = rmfield(AllMeasurements,'parameters');
    end
    handles.Measurements = AllMeasurements;
    set(handles.text2A,'String',num2str(size(handles.Measurements.CellData,1)));
    set(handles.text2C,'String',num2str(length(unique(handles.Measurements.CellData(:,1)))));
    set(handles.text2E,'String',num2str(max(handles.Measurements.CellData(:,4))));
    set(handles.text3B,'String','--')
    % Make experiment name (for graph titles) from folder structure
    folder = get(handles.edit1A,'String');
    ind1 = regexp(folder,'20..-..-');
    ind2 = strfind(folder,filesep);
    ind2 = ind2(end);
    if ~isempty(ind1)
        date1 = folder(ind1:ind2-1);
        condition1 = folder(ind2+1:end);
        handles.Export.Name = [condition1,' - ', date1];
        handles.Export.Name(strfind(handles.Export.Name,'_')) = ' ';
        disp('- - - - - - - - - -')
        disp(['Loading measurements: "',handles.Export.Name,'"'])
        num1 = sum( (handles.Measurements.CellData(:,3)<13)&(handles.Measurements.CellData(:,4)>11));
        num2 = 0.1*round(10*num1/length(unique(handles.Measurements.CellData(:,1))));
        disp(['-> ',num2str(num1),' cells at 12th timepoint (',num2str(num2),' cells per XY)'])
        disp('- - - - - - - - - -')

    else
        handles.Export.Name = '';
    end
    
    
else
    disp('Invalid Measurements file')
    set(handles.text2A,'String','--');
    set(handles.text2C,'String','--');
    set(handles.text2E,'String','--');
    set(handles.text3B,'String','--')
end

guidata(handles.figure1,handles)
% ========================================================================================

% - - - - - - - - - - - - - - - - - - - - UIPANEL 3: CELL FILTERING  - - - - - - - - - - - - - - - - - - - -

function checkbox3A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX3A: filter by non-edge cells
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.FilterFlags(1) = 1;
else
    handles.Options.FilterFlags(1) = 0;
end
guidata(handles.figure1,handles);
% ========================================================================================

function checkbox3B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX3B:filter by parents (cell division) only
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.FilterFlags(2) = 1;
else
    handles.Options.FilterFlags(2) = 0;
end
guidata(handles.figure1,handles);
% ========================================================================================

function checkbox3C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX3C:Filter by dell death
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.FilterFlags(3) = 1;
else
    handles.Options.FilterFlags(3) = 0;
end
guidata(handles.figure1,handles);
% ========================================================================================

function checkbox3D_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX3D:filter by minimum cell lifetime (in frames)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.FilterFlags(4) = 1;
else
    handles.Options.FilterFlags(4) = 0;
end
guidata(handles.figure1,handles);
% ========================================================================================

function edit3A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT3A: set minimum lifetime for checkbox3D
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.MinLifetime = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);

% ========================================================================================

function pushbutton3A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON3A: filter cells by selected criteria, update text
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

[DropCells] = filterCells(handles.Measurements.CellData, handles.Options.FilterFlags, handles.Options.MinLifetime);
handles.DropCells = DropCells;    
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1,handles);
set(handles.text3B,'String',num2str(size(handles.Measurements.CellData,1)-sum(handles.DropCells)));
loadmeasurefields(handles)
% ========================================================================================

function loadmeasurefields(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Load measurements when cells are filtered (subfunction)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

measurementNames = fieldnames(handles.Measurements);

% Drop out any cell matricies (we're not equipped to graph them)
dropList = [];
for i = 1:length(measurementNames)
    if iscell(handles.Measurements.(measurementNames{i}))
        dropList = cat(1,dropList,i);
    end
end
measurementNames(dropList) = [];

measurementNames = measurementNames(~strcmp(measurementNames,'CellData'));
% Update MeasureField1
set(handles.listbox4A,'String',measurementNames,...
 'Value',1)
handles.Export.MeasureField1 = measurementNames{1};
% Update MeasureField2
set(handles.listbox5A,'String',measurementNames,...
 'Value',1)
handles.Export.MeasureField2 = measurementNames{1};
guidata(handles.figure1,handles);
% ========================================================================================



% - - - - - - - - - - - - - - - - - - - - UIPANEL 4-5: MEASUREMENT SELECTION  - - - - - - - - - - - - - - - - - - - -

function listbox4A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX4A: save selected measurement to 'handles'
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

index_selected = get(handles.listbox4A,'Value');
measure_list = get(handles.listbox4A,'String');
% Update MeasureField1 to new selected item
if ~strcmp(handles.Export.MeasureField1,measure_list{index_selected})
    handles.Export.MeasureField1 = measure_list{index_selected};
    handles.Flags.Modeled = 0; % Reset "modeled" flag
    guidata(handles.figure1,handles);
end
% ========================================================================================

function listbox5A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LISTBOX5A: save selected measurement to 'handles'
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

index_selected = get(handles.listbox5A,'Value');
measure_list = get(handles.listbox5A,'String');
% Update MeasureField2 to new selected item
if ~strcmp(handles.Export.MeasureField2,measure_list{index_selected})
    handles.Export.MeasureField2 = measure_list{index_selected};
    handles.Flags.Modeled = 0; % Reset "modeled" flag
    guidata(handles.figure1,handles);
end
%========================================================================================




% - - - - - - - - - - - - - - - - - - - - UIPANEL 6: DATA SMOOTHING  - - - - - - - - - - - - - - - - - - - -

function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BUTTONGROUP: record smoothing option to 'handles'
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if ~strcmp(handles.Options.Smoothing,get(eventdata.NewValue,'String'))
    handles.Options.Smoothing = get(eventdata.NewValue,'String');
    handles.Flags.Modeled = 0; % Reset "modeled" flag
    guidata(handles.figure1, handles)
end
% ========================================================================================

function edit6A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT6A: set sliding window size for data smoothing
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
val = floor(str2double(get(hObject,'String')));
set(hObject,'String',num2str(val))
if (handles.Options.SmoothingWindow~= val)
    handles.Options.SmoothingWindow = val;
    handles.Flags.Modeled = 0; % Reset "modeled" flag
    guidata(handles.figure1,handles);
end
% ========================================================================================



% - - - - - - - - - - - - - - - - - - - - UIPANEL 7: DATA OPTIONS  - - - - - - - - - - - - - - - - - - - -

function edit7A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7A: set pixel-to-micron calibration
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.PixelConversion = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================

function edit7B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7B: set frames-per-hour of image data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.FramesPerHour = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================

function edit7C_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7C: manually set graphed minimum for Measurement1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.GraphMin1 = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================

function edit7D_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7D: manually set graphed maximum for Measurement2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.GraphMin2 = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================



function edit7E_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7C: manually set graphed maximum for Measurement1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.GraphMax1 = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================



function edit7F_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7D: manually set graphed maximum for Measurement2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.GraphMax2 = str2double(get(hObject,'String'));
guidata(handles.figure1,handles);
% ========================================================================================


function checkbox7A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX7A: log compression of Measurement1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.LogCompress1 = 1;
else
    handles.Options.LogCompress1 = 0;
end
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1,handles);
% ========================================================================================

function checkbox7B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX7B: log compression of Measurement2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.LogCompress2 = 1;
else
    handles.Options.LogCompress2 = 0;
end
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1,handles);
% ========================================================================================



function edit7G_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7G: set number of bins for histogram 1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.HistogramBins1 = round(str2double(get(hObject,'String')));
set(hObject,'String', num2str(handles.Options.HistogramBins1))
guidata(handles.figure1,handles);
% ========================================================================================


function edit7H_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7H: set number of bins for histogram 2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.HistogramBins2 = round(str2double(get(hObject,'String')));
set(hObject,'String', num2str(handles.Options.HistogramBins2))
guidata(handles.figure1,handles);
% ========================================================================================

function edit7I_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7I: set histogram1 maximum (Options.HistogramMax1)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.HistogramMax1 = str2double(get(hObject,'String'));
set(hObject,'String', num2str(handles.Options.HistogramMax1))
guidata(handles.figure1,handles);
% ========================================================================================

function edit7J_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7J: set histogram2 maximum (Options.HistogramMax2)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.HistogramMax2 = str2double(get(hObject,'String'));
set(hObject,'String', num2str(handles.Options.HistogramMax2))
guidata(handles.figure1,handles);
% ========================================================================================


function edit7K_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7K: set time-axis start time hours (Options.StartTime)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

timeHr = str2double(get(handles.edit7K,'String'));
timeMin = str2double(get(handles.edit7K_2,'String'));

handles.Options.StartTime = timeHr+timeMin/60;
set(handles.edit7K,'String',numseq(fix(handles.Options.StartTime),2))
set(handles.edit7K_2,'String',numseq(round((handles.Options.StartTime-fix(handles.Options.StartTime))*60),2))
guidata(handles.figure1,handles);
% ========================================================================================

function edit7K_2_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT7K_2: set time-axis start time minutes (Options.StartTime)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

timeHr = str2double(get(handles.edit7K,'String'));
timeMin = str2double(get(handles.edit7K_2,'String'));

handles.Options.StartTime = timeHr+timeMin/60;
set(handles.edit7K,'String',numseq(fix(handles.Options.StartTime),2))
set(handles.edit7K_2,'String',numseq(round((handles.Options.StartTime-fix(handles.Options.StartTime))*60),2))
guidata(handles.figure1,handles);
% ========================================================================================


% - - - - - - - - - - - - - - - - - - - - UIPANEL 8: Grouping/Statistical Modeling  - - - - - - - - - - - - - - - - - -

function uipanel8_SelectionChangeFcn(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BUTTONGROUP8: select new statistical modeling type
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.Options.Grouping = get(eventdata.NewValue,'String');
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1, handles)
% ========================================================================================

function edit8A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT8A: set number of subpopulations used in Gaussian modeling (Options.Subpopulations)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.Options.Subpopulations = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(handles.Options.Subpopulations))
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1,handles);
% ========================================================================================

function checkbox8A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKBOX8A: if selected, group using measurement 2 (Options.Measure2Grouping)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (get(hObject,'Value')==get(hObject,'Max'))
    handles.Options.Measure2Grouping = 1;
else
    handles.Options.Measure2Grouping = 0;
end
handles.Flags.Modeled = 0; % Reset "modeled" flag
guidata(handles.figure1,handles);
% ========================================================================================


% - - - - - - - - - - - - - - - - - - - - UIPANEL 9: GRAPHING AND EXPORT  - - - - - - - - - - - - - - - - - - - -

function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BUTTONGROUP9: record selected graph type to 'handles'
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.GraphType = get(eventdata.NewValue,'String');
guidata(handles.figure1, handles)
% ========================================================================================

function pushbutton9A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON9A: call 'PlotMeasurement' to graph with selected options in figure 1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if ~handles.Flags.Modeled
    handles = processMeasurements(handles);
end
handles.Flags.Modeled = 1; % Set "modeled" flag  

handles = setGraphOptions(handles);
guidata(handles.figure1, handles)

% If hierarchial clustering is on, create additional dendrogram
if strcmp(get(handles.radiobutton8D,'String'),handles.Options.Grouping)
    handles.AltFigure = figure;
    set(handles.AltFigure,'DockControls','off','Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')
    C = cell(size(handles.Export.Measurement1,1),1);
    C(:) = {''};
    figure(handles.AltFigure),[h, ~,perm] = dendrogram(handles.Export.Linkage,0,'orientation','left','labels',C);
    dendrogram_axes = get(h(1),'Parent');
    set(dendrogram_axes,'LooseInset',get(dendrogram_axes,'TightInset')+ [0 0.03 0.02 0.01])
    set(h,'Color',[0.3 0.3 0.3])
end


% Parse graph type and plot
export = handles.Export;
handles.PlotFigure = figure;% Create new figure, set properties
set(0,'defaultAxesFontName', 'arial')
set(handles.PlotFigure,'Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')
switch handles.Options.GraphType
    case get(handles.radiobutton9A,'String') %Line plot
        linePlot(export.Measurement1,export.CellData,export.Measurement1Info,export.GroupingVector, handles.PlotFigure)
    case get(handles.radiobutton9B,'String') % Stack colormap
        colormapStack(export.Measurement1,export.CellData,export.Measurement1Info,handles.PlotFigure);
    case get(handles.radiobutton9C,'String') % Histogram series
        delete(get(handles.PlotFigure,'Children'))
        histogramSeries(handles);
    case get(handles.radiobutton9D,'String') % Scatterplot
        delete(get(handles.PlotFigure,'Children'))
        scatterPlot(handles);
end
% Assignin variables into workspace.
assignin('base', 'handles',handles)
assignin('base', 'export', handles.Export)
% ========================================================================================

% --- Executes on button press in pushbutton9B.
function pushbutton9B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON9B: call 'ExportData' to save filtered/smoothed data for further manipulation
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
assignin('base', 'handles',handles)
assignin('base', 'export', handles.Export)
% ========================================================================================


% - - - - - - - - - - - - - - - - - - - - UIPANEL 10: CELL VISUALIZATION  - - - - - - - - - - - - - - - - - - - -

function edit10A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT10A: set xy position to be visualized outlined on images
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

handles.Options.VisualizeXY = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(handles.Options.VisualizeXY))
guidata(handles.figure1, handles)


function edit10B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT10B: set number of cell to be visualized outline on images
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.Options.VisualizeCell = round(str2double(get(hObject,'String')));
set(hObject,'String',num2str(handles.Options.VisualizeCell))
guidata(handles.figure1, handles)

function pushbutton10A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON10A: create new figure and call visualization function
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
visualizeCell(handles)

function pushbutton10B_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON10A: create new figure and call visualization function
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ~exist(handles.locations.data,'dir')
    specifyLocations;
    uiwait;
    load([handles.home_folder, 'locations.mat'],'-mat')
    handles.locations = locations;
    guidata(handles.figure1,handles);
end
checkDynamics(handles.file_loc);


% - - - - - - - - - - - - - - - - - - - - UIPANEL 11: XY VISUALIZATION  - - - - - - - - - - - - - - - - - - - -
function pushbutton11A_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON11A: call XY visualization function
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ~exist(handles.locations.data,'dir')
    specifyLocations;
    uiwait;
    load([handles.home_folder, 'locations.mat'],'-mat')
    handles.locations = locations;
    guidata(handles.figure1,handles);
end
showTiles([handles.locations.scope, handles.trackparameters.ImagePath]);

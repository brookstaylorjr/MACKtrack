function varargout = specifyLocations(varargin)
% SPECIFYLOCATIONS MATLAB code for specifyLocations.fig
%      SPECIFYLOCATIONS, by itself, creates a new SPECIFYLOCATIONS or raises the existing
%      singleton*.
%
%      H = SPECIFYLOCATIONS returns the handle to a new SPECIFYLOCATIONS or the handle to
%      the existing singleton*.
%
%      SPECIFYLOCATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECIFYLOCATIONS.M with the given input arguments.
%
%      SPECIFYLOCATIONS('Property','Value',...) creates a new SPECIFYLOCATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before specifyLocations_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to specifyLocations_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help specifyLocations

% Last Modified by GUIDE v2.5 13-Apr-2015 14:33:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @specifyLocations_OpeningFcn, ...
                   'gui_OutputFcn',  @specifyLocations_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


function specifyLocations_OpeningFcn(hObject, eventdata, handles, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to secifyLocations (see VARARGIN)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Choose default command line output for specifyLocations
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set colors
handles.blue = [0.051 0.329 0.671];
handles.orange = [207 79 51]/255;
handles.gray = [0.816 0.812 0.722];


% Load existing tracking.mat structure and populate edit fields
p = mfilename('fullpath');
slash_idx = strfind(p,filesep);
handles.save_loc = [p(1:slash_idx(end-1)), 'locations.mat'];
load(handles.save_loc,'-mat')
if ~strcmp(locations.scope(end),filesep)
    locations.scope = [locations.scope,filesep];
end
if ~strcmp(locations.data(end),filesep)
    locations.data = [locations.data,filesep];
end
handles.locations = locations;

set(handles.edit1,'String',handles.locations.scope,'ForegroundColor',handles.gray);
set(handles.edit2,'String',handles.locations.data,'ForegroundColor',handles.gray);
set(handles.edit3,'String',handles.locations.spreadsheet,'ForegroundColor',handles.blue);
handles.Locked1 = 1;
handles.Locked2 = 1;
handles.Locked3 = 0;
if exist(handles.locations.scope,'dir') 
    handles.Locked1 = 0;
    set(handles.edit1,'ForegroundColor',handles.blue)
end
if  exist(handles.locations.data,'dir')
    handles.Locked2 = 0;
    set(handles.edit2,'ForegroundColor',handles.blue)
end

if ~handles.Locked1 && ~handles.Locked2
    set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
end

guidata(handles.figure1,handles)


function varargout = specifyLocations_OutputFcn(hObject, eventdata, handles) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set default command line output from handles structure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT1: set location of mount point for device that holds scope images (this way,
% multiple computers can access the same images without an error)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newfolder = get(hObject,'String');
    if exist(newfolder,'dir')
        if ~strcmp(newfolder(end),filesep)
            newfolder = [newfolder,filesep];
        end
        set(handles.edit1,'String',newfolder,'ForegroundColor',handles.blue)
        handles.locations.scope = newfolder;
        handles.Locked1 = 0;
        if ~handles.Locked2
            set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
        end
        guidata(handles.figure1,handles)
    else
        handles.Locked1 = 1;
        set(handles.edit1,'ForegroundColor',handles.gray)
        set(handles.pushbutton_end,'ForegroundColor',handles.gray)
        guidata(handles.figure1,handles)
    end
catch ME
    set(hObject,'String','err')
    handles.Locked = 1;
    set(handles.edit1,'ForegroundColor',handles.gray)
    set(handles.pushbutton_end,'ForegroundColor',handles.gray)
    guidata(handles.figure1,handles)
    rethrow(ME)
end
% ========================================================================================


function pushbutton1_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1: browse for folder (scope data)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
PathName = uigetdir();
if ~strcmp(PathName(end),filesep)
	PathName = [PathName,filesep];
end
if (PathName~=0)
    set(handles.edit1,'String',PathName)
    set(handles.edit1,'ForegroundColor',handles.blue)
    handles.locations.scope = PathName;
    % Unlock save button
    handles.Locked1 = 0;
    if ~handles.Locked2
        set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
    end
    guidata(handles.figure1,handles)
end
% ========================================================================================

function edit2_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT2: set location of mount point for device that holds data (this way,
% multiple computers can access the same datasets without an error)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newfolder = get(hObject,'String');
    if exist(newfolder,'dir')
        if ~strcmp(newfolder(end),filesep)
            newfolder = [newfolder,filesep];
        end
        set(handles.edit2,'String',newfolder, 'ForegroundColor',handles.blue)
        handles.locations.data = newfolder;
        handles.Locked2 = 0;
        if ~handles.Locked1
            set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
        end
        guidata(handles.figure1,handles)
    else
        handles.Locked2 = 1;
        set(handles.edit2,'ForegroundColor',handles.gray)
        set(handles.pushbutton_end,'ForegroundColor',handles.gray)
        guidata(handles.figure1,handles)
    end
catch ME
    set(hObject,'String','err')
    handles.Locked2 = 1;
    set(handles.edit1,'ForegroundColor',handles.gray)
    set(handles.pushbutton_end,'ForegroundColor',handles.gray)
    guidata(handles.figure1,handles)
    rethrow(ME)
end
% ========================================================================================




function pushbutton2_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON1: browse for folder
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
PathName = uigetdir();
if ~strcmp(PathName(end),filesep)
	PathName = [PathName,filesep];
end
if (PathName~=0)
    set(handles.edit2,'String',PathName)
    set(handles.edit2,'ForegroundColor',handles.blue)
    handles.locations.data = PathName;
    handles.Locked2 = 0;
    % Unlock save button
    if ~handles.Locked1
        set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
    end
    guidata(handles.figure1,handles)
end
% ========================================================================================


function edit3_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDIT3: set URL of Google Spreadsheet that holds 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try
    newaddr = get(hObject,'String');
    handles.locations.spreadsheet = newaddr;
    guidata(handles.figure1,handles)
    set(handles.edit3,'ForegroundColor',handles.blue)
    set(handles.pushbutton_end,'ForegroundColor',[0 0 0])
    handles.Locked3 = 0;

catch ME
    set(hObject,'String','err')
    handles.Locked3 = 1;
    set(handles.edit1,'ForegroundColor',handles.gray)
    set(handles.pushbutton_end,'ForegroundColor',handles.gray)
    guidata(handles.figure1,handles)
    rethrow(ME)
end
% ========================================================================================


function pushbutton_end_Callback(hObject, eventdata, handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PUSHBUTTON_END: save to 'locations.mat'
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ~handles.Locked1 && ~handles.Locked2 && ~handles.Locked3
    locations = handles.locations;
    save(handles.save_loc,'locations')
    close(handles.figure1);
end
% ========================================================================================

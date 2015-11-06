function varargout = checkDynamics(varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = checkDynamics(id)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKDYNAMICS creates a figure showing tracked/outlined cells that corrrespond to a 
% dynamic measurement, graphed below. Measurements are first processed by specific 
% functions called (e.g.) see_nfkb or see_nfkb_dim.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%  MATLAB code generated for checkDynamics.fig:
%      checkDynamics, by itself, creates a new checkDynamics or raises the existing
%      singleton*.
%
%      H = checkDynamics returns the handle to a new checkDynamics or the handle to
%      the existing singleton*.
%
%      checkDynamics('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in checkDynamics.M with the given input arguments.
%
%      checkDynamics('Property','Value',...) creates a new checkDynamics or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before checkDynamics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to checkDynamics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help checkDynamics

% Last Modified by GUIDE v2.5 17-Feb-2015 11:57:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkDynamics_OpeningFcn, ...
                   'gui_OutputFcn',  @checkDynamics_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if (nargin>1) && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before checkDynamics is made visible.
function checkDynamics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to checkDynamics (see VARARGIN)

% Choose default command line output for checkDynamics
handles.output = hObject;

% Initialize dropdown menu with visualization scripts; chose translocation by default
home_folder = mfilename('fullpath');
home_folder = home_folder(1:max(strfind(home_folder,filesep)));
fcnlist = dir(home_folder);
viz_fcns = cell(0);
for i = 1:length(fcnlist)
    if strcmp(fcnlist(i).name(1:min(end,4)),'see_')
        viz_fcns = cat(1,viz_fcns,fcnlist(i).name(1:end-2));
    end
end
set(handles.popupmenu2,'String',viz_fcns,'Value',1)
handles.DataFcn = viz_fcns(1);


% Get directory locations - locations.mat is one directory up.
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
if ~exist(locations.scope,'dir') || ~exist(locations.data,'dir')
    specifyLocations;
    uiwait;
    load([home_folder(1:slash_idx(end)), 'locations.mat'],'-mat')
end
handles.locations = locations;

% Load visualization data into handles
handles.id = varargin{1};
handles = load_vizdata(handles);
guidata(handles.figure1, handles)


hZoom = zoom(gcf);
set(hZoom,'ActionPostCallback',{@customZoom,handles});


% Set output function
function varargout = checkDynamics_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function handles_out = load_vizdata(handles)
fcn_names = get(handles.popupmenu2,'String');
start_val = get(handles.popupmenu2,'Value');
max_val = length(fcn_names);
id = handles.id;
% Read in function from dropdown menu, and process data appropriately
disp('loading data...')
flag = 0;
while flag == 0
    try
        fcn_name = fcn_names{start_val};
        [graph, info] = eval([fcn_name,'(id,0,0)']);
        set(handles.popupmenu2,'Value',start_val);
        flag = 1;
    catch ME
        start_val = start_val + 1;
        disp(['Skipped function "' fcn_name,'"'])
        disp('- - - - - - (begin error trace) - - - - - - - - - - - - -')
        disp(getReport(ME,'extended'));
        disp('- - - - - - (end error trace) - - - - - - - - - - - - - -')
        if start_val>max_val
            error('No valid measurement funtions found.')
        end
        flag = 0;
    end
end

handles.times = info.parameters.TimeRange;
handles.xys = info.parameters.XYRange;
handles.celldata = graph.celldata;
handles.var = graph.var;
handles.mu = nanmean(handles.var);
handles.sigma = nanstd(handles.var);
handles.t = graph.t;
handles.shift = graph.shift;
handles.ylim = info.graph_limits;
handles.parameters = info.parameters;
handles.module = info.Module;
if ~isempty(strfind(info.savename,filesep))
    handles.mask_dir = info.savename(1:max(strfind(info.savename,filesep)));
else
    handles.mask_dir = [pwd,filesep];
end

% Initialize slider + popup values
handles.xys =  info.parameters.XYRange;
set(handles.slider1, 'Min',1, 'Max',length(handles.xys),'SliderStep',[1 4]/length(handles.xys),'Value',1);
handles.xy = handles.xys(1);
set(handles.text5,'String',['xy ',num2str(handles.xy)]);
v1 = 1; v2 = length(handles.t); 
set(handles.slider2, 'Min',v1, 'Max',v2,'SliderStep',[1/(v2-v1) 4/(v2-v1)],'Value',v1);
handles.time = v1;
% Initialize image+graph, update handles
handles = newXY(handles);
handles_out = handles;





% GUI callbacks: update values, call functions
function slider1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderval = round(get(hObject,'Value'));
handles.xy = handles.xys(sliderval);
set(hObject,'Value',sliderval);
set(handles.text5,'String',['xy ',num2str(handles.xy)]);
handles = newXY(handles);
guidata(handles.figure1, handles);


function slider2_Callback(hObject, eventdata, handles)
handles.time = round(get(hObject,'Value'));
set(hObject,'Value',handles.time);
handles = newTime(handles);
guidata(handles.figure1, handles);

function popupmenu1_Callback(hObject, eventdata, handles)
cellno = handles.cell_list(get(handles.popupmenu1,'Value'));
if cellno ~= handles.cell
    handles.cell = cellno;
    handles = newCell(handles);
end
guidata(handles.figure1, handles);


function popupmenu2_Callback(hObject, eventdata, handles)
% Change function in use - on selection change, recalculate and reassign measurements into handles.
% Load visualization data into handles
handles = load_vizdata(handles);
guidata(handles.figure1, handles);


function customZoom(fig,evd,handles) %#ok<INUSL>
% Modify zoom behavior to make sure that axes are reset to full image on double click (either zoom in or zoom out)
if strcmp(get(fig,'SelectionType'),'open')
    xlim(handles.axes1,[1 size(handles.nuc_label,2)])
    ylim(handles.axes1,[1 size(handles.nuc_label,1)])
end

% Update information in axes1 and axes2 as needed
function handles = newXY(handles)
% Make the list of cells in the current xy, update popup menu
handles.cell_list = handles.celldata(handles.celldata(:,1)==handles.xy,2);
set(handles.popupmenu1,'String',cellstr(num2str(handles.cell_list)),'Value',1)
handles.cell = handles.cell_list(get(handles.popupmenu1,'Value'));


% Update and draw image
handles = loadImage(handles);
drawImage(handles)
drawGraph(handles)

function handles = newTime(handles)
% Update and draw image
handles = loadImage(handles);
drawImage(handles);
drawGraph(handles)


function handles = newCell(handles)
% Update image (mask only)
drawImage(handles);
% Update graph
drawGraph(handles)


function handles = loadImage(handles)
% Get new image/masks
i = handles.xy;
xy_idx = find(handles.xys==handles.xy,1,'first');
j = handles.time-handles.shift(xy_idx);
img_path = [handles.locations.scope, filesep, handles.parameters.ImagePath,filesep,...
    eval(handles.parameters.(handles.module).ImageExpr)];
% Get bit depth of image, load, and rotate (if necessary)
if ~isfield(handles.parameters,'BitDepth')
    imfo = imfinfo(img_path);
    handles.parameters.BitDepth = imfo.BitDepth;
end
measure_img = checkread(img_path,handles.parameters.BitDepth);
if size(measure_img,1)>size(measure_img,2)
    measure_img = imrotate(measure_img,90);
    flip_flag = 1;
end
% Load masks, rotate if necessary
load([handles.mask_dir,filesep,'xy',num2str(i),filesep,...
    'NuclearLabels',filesep,'NuclearLabel-',numseq(j,4),'.mat'])
load([handles.mask_dir,'xy',num2str(i),filesep,...
    'CellLabels',filesep,'CellLabel-',numseq(j,4),'.mat'])
if exist('flip_flag','var')
    handles.nuc_label = imrotate(NuclearLabel,90);
    handles.cell_label = imrotate(CellLabel,90);
else
    handles.nuc_label = NuclearLabel;
    handles.cell_label = CellLabel;
end


% Saturate image according to first image
if 1%~isfield(handles,'imgmax')
handles.imgmax = prctile(measure_img(:),99.9);
handles.imgmin = prctile(measure_img(:),5);
end
measure_img(measure_img>handles.imgmax) = handles.imgmax;
measure_img(measure_img<handles.imgmin) = handles.imgmin;
handles.img = uint8((measure_img-handles.imgmin)/(handles.imgmax-handles.imgmin)*255);
% Get centroids and labels
props = regionprops(handles.nuc_label,'Centroid');
tmpcell = struct2cell(props);
tmpmat = cell2mat(tmpcell(1,:));
handles.centroids = [tmpmat(1:2:end)', tmpmat(2:2:end)'];
handles.centroids(isnan(handles.centroids(:,1)),:) = [];
handles.text_labels = unique(NuclearLabel(NuclearLabel>0));
handles.centroids(~ismember(handles.text_labels,handles.cell_list),:)=[];
handles.text_labels = cellstr(num2str(handles.text_labels(ismember(handles.text_labels,handles.cell_list))));



function drawImage(handles)
% Color masks
if ~isempty(get(handles.axes1,'Children'))
    xlim = get(handles.axes1,'XLim');
    ylim = get(handles.axes1,'YLim');
    setflag = 1;
else
    setflag = 0;  
end
a = 0.24;
img_a = handles.img; img_b = handles.img; img_c = handles.img;
mask1 = ismember(handles.cell_label,handles.cell_list);
borders = handles.cell_label;
borders(~mask1) = 0;
borders = (imdilate(borders,true(5)) - borders)>0;
mask1 = mask1&~ismember(handles.nuc_label,handles.cell_list);
img_a(mask1) = handles.img(mask1)*(1-a) + 45*a;
img_b(mask1) = handles.img(mask1)*(1-a) + 103*a;
img_c(mask1) = handles.img(mask1)*(1-a) + 207*a;

a2 = 0.6;
img_a(borders) = handles.img(borders)*(1-a2) + 45*a2;
img_b(borders) = handles.img(borders)*(1-a2) + 103*a2;
img_c(borders) = handles.img(borders)*(1-a2) + 207*a2;

% Color current cell
mask2 = (handles.cell_label==handles.cell) &~(handles.nuc_label==handles.cell);
img_a(mask2) = handles.img(mask2)*(1-a) + 232*a;
img_b(mask2) = handles.img(mask2)*(1-a) + 119*a;
img_c(mask2) = handles.img(mask2)*(1-a) + 44*a;
borders2 = imdilate(mask2,true(5))&~mask2;
img_a(borders2) = handles.img(borders2)*(1-a2) + 232*a2;
img_b(borders2) = handles.img(borders2)*(1-a2) + 119*a2;
img_c(borders2) = handles.img(borders2)*(1-a2) + 44*a2;
imshow(cat(3,img_a,img_b,img_c),'Parent',handles.axes1)
set(handles.axes1,'XTick',[],'YTick',[])
if setflag
    set(handles.axes1,'xlim',xlim,'ylim',ylim,'XTick',[],'YTick',[])
else
    axis image
end
text(handles.centroids(:,1),handles.centroids(:,2),handles.text_labels,'FontSize',14,...
    'Color',[219 85 31]/255,'Parent',handles.axes1,'HorizontalAlignment','center')


function drawGraph(handles)
mu_xy = nanmean(handles.var(handles.celldata(:,1)==handles.xy,:));
sigma_xy = nanstd(handles.var(handles.celldata(:,1)==handles.xy,:));

fill([handles.t,handles.t(end:-1:1)],[handles.mu+handles.sigma,handles.mu(end:-1:1)-handles.sigma(end:-1:1)],...
    [.8 .8 .8],'Parent',handles.axes2,'FaceAlpha',0.5)
idx = handles.celldata(:,1)==handles.xy&handles.celldata(:,2)==handles.cell;
hold(handles.axes2,'on')
fill([handles.t,handles.t(end:-1:1)],[mu_xy+sigma_xy,mu_xy(end:-1:1)-sigma_xy(end:-1:1)],[68 129 227]/255,...
    'FaceAlpha',0.25,'Parent',handles.axes2)
plot(handles.axes2,handles.t, handles.var(idx,:),'Color',[68 129 227]/255,'LineWidth',3)
hold(handles.axes2,'off')
axis(handles.axes2,[min(handles.t),max(handles.t),handles.ylim])
handles.line1 = line([handles.t(handles.time) handles.t(handles.time)],handles.ylim,...
    'Color',[0 0 0],'Parent',handles.axes2);
handles.line2 = line([min(handles.t) max(handles.t)],[handles.var(idx,handles.time) handles.var(idx,handles.time)],...
    'Color',[0 0 0],'Parent',handles.axes2);
set(handles.axes2,'Layer','Top')
y_lim = get(handles.axes2,'YLim');
text(max(handles.t)-range(handles.t)*0.12,max(y_lim)-range(y_lim)*0.07,...
    ['nuclear NF\kappaB: ',num2str(.01*round(handles.var(idx,handles.time)*100))],'Parent',handles.axes2)

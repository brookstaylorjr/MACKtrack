function varargout = checkDynamics(varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = checkDynamics(varargin)
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

% Last Modified by GUIDE v2.5 17-Feb-2015 11:57:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @checkDynamics_OpeningFcn, ...
                   'gui_OutputFcn',  @checkDynamics_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
if nargin<1
    error('Input must include an experiment ID/file location')
end
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


% 1) Get list of visualization functions (named "see_...")
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

% 2) Get system-specific locations
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
if ~exist(locations.scope,'dir') || ~exist(locations.data,'dir')
    setLocations;
    uiwait;
    load([home_folder(1:slash_idx(end)), 'locations.mat'],'-mat')
end
handles.locations = locations;

% 3) Load specified measurements file - get source image directory and output label matricies
[~, info, handles.AllMeasurements] = loadID(varargin{1});
handles.ImageDir = namecheck(info.ImageDirectory);
if ~exist(handles.ImageDir,'dir')
    error(['Could not find specified image directory, ', handles.ImageDir])
end
handles.OutputDir = 'notadir';
if ischar(varargin{1})
    fullname = varargin{1};
    slash_idx = strfind(fullname,filesep);
    if exist(fullname(1:slash_idx(end)),'NuclearLabels','dir')
        handles.OutputDir = fullname(1:slash_idx(end));
    end
end
if ~exist(handles.OutputDir,'dir')
    handles.OutputDir = namecheck([locations.data,filesep,handles.AllMeasurements.parameters.SaveDirectory]);
    if ~exist(handles.OutputDir,'dir')
        error(['Could not find specified output directory, ', handles.OutputDir])
    end
end


% Load visualization data into handles
set(handles.popupmenu2,'Value',1)
handles.FcnVal = [];
handles = load_vizdata(handles);
guidata(handles.figure1, handles)

hZoom = zoom(gcf);
set(hZoom,'ActionPostCallback',{@customZoom,handles});


% Set output function
function varargout = checkDynamics_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function handles_out = load_vizdata(handles)
fcn_names = get(handles.popupmenu2,'String');
newval = get(handles.popupmenu2,'Value');
max_val = length(fcn_names);
% Read in function from dropdown menu, and process data appropriately
    flag = 0;
    while flag == 0
        try
            fcn_name = fcn_names{newval};
            disp('- - - - ')
            disp(['Trying to use display function ', upper(fcn_name)])
            [graph, info] = eval([fcn_name,'(handles.AllMeasurements)']);
            set(handles.popupmenu2,'Value',newval);
            handles.FcnVal = newval;
            flag = 1;
        catch ME
            if isempty(handles.FcnVal)
                newval = newval + 1;
                disp(['Skipping - hit error:'])
                disp(getReport(ME));
                if newval>max_val
                    error('No valid measurement funtions found.')
                end
                flag = 0;
            else
                set(handles.popupmenu2,'Value', handles.FcnVal)
                disp('Invalid - hit error:')
                disp(getReport(ME,'basic'));
                fcn_name = fcn_names{handles.FcnVal};
                disp(['Resetting to ', upper(fcn_name)])
                [graph, info] = eval([fcn_name,'(handles.AllMeasurements)']);
                flag = 1;
            end
        end
    end


% Get primary fields- timepoints & data
handles.times = info.parameters.TimeRange;
handles.xys = unique(graph.celldata(:,1));
handles.celldata = graph.celldata;
handles.var = graph.var;
handles.t = graph.t;
handles.parameters = info.parameters;
handles.ImageExpr = info.ImageExpr;

% Add other fields - division matrix, 2nd measurement, etc. 
if isfield(graph,'lineage'); 
    handles.lineage = graph.lineage; 
else
    handles.lineage = handles.celldata(:,2);
end
if isfield(graph,'var2'); handles.var2 = graph.var2; end
if isfield(graph,'shift'); 
    handles.shift = graph.shift; 
else
    handles.shift = zeros(length(unique(handles.celldata(:,1))),1);
end
if isfield(info,'graph_limits'); 
    handles.ylim = info.graph_limits; 
else
    handles.ylim = prctile(handles.var(:),[2 98]);
end

% Make calculated fields - average across all cells per timept
handles.mu = nanmean(handles.var);
handles.sigma = nanstd(handles.var);


% Initialize slider + popup values
set(handles.slider1, 'Min',1, 'Max',length(handles.xys),'SliderStep',[1 4]/length(handles.xys),'Value',1);
handles.xy = handles.xys(1);
set(handles.text5,'String',['xy ',num2str(handles.xy)]);
v1 = 1; v2 = length(handles.t); 
set(handles.slider2, 'Min',v1, 'Max',v2,'SliderStep',[1/(v2-v1) 4/(v2-v1)],'Value',v1);
handles.time = v1;
% Initialize image+graph, update handles
handles = newXY(handles);
handles_out = handles;
guidata(handles.figure1, handles)


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
if cellno ~= handles.cell_idx
    handles.cell_idx = cellno;
    guidata(handles.figure1, handles);
    handles = newCell(handles);
end


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


% - - - - - GUI UPDATE functions: - - - - - - - - - - - - - - - - - - - - - - - - - - - 
function handles = newXY(handles)
% Get subset + list of cells in the current xy, update popup menu
handles.xy_subset = find(handles.celldata(:,1)==handles.xy);
handles.cell_list = cell(size(handles.xy_subset));
for i = 1:length(handles.cell_list)
    tmp_row = unique(handles.lineage(handles.xy_subset(i),:));
    tmp_row(isnan(tmp_row)) = [];
    tmp_row = num2str(tmp_row);
    tmp_row(strfind(tmp_row,'  ')) = '-';
    tmp_row(strfind(tmp_row,' ')) = '>';
    handles.cell_list{i} = tmp_row;
end
set(handles.popupmenu1,'String',cellstr(handles.cell_list),'Value',1)
handles.cell_idx = get(handles.popupmenu1,'Value');
handles.cells_all = unique(handles.lineage(handles.xy_subset,:));
handles.cells_all(isnan(handles.cells_all)) = [];

% Create basic metrics for subset
handles.mu_subset = nanmean(handles.var(handles.xy_subset,:));
handles.sigma_subset = nanstd(handles.var(handles.xy_subset,:));

% New XY: load image/labelmats, mask/display cells, graph trajectory
handles = loadImage(handles);
drawImage(handles)
drawGraph(handles)


function handles = newTime(handles)
% New time: load image/labelmats, mask/display cells, graph trajectory
handles = loadImage(handles);
drawImage(handles);
drawGraph(handles)


function handles = newCell(handles)
% New time: mask/display cells, graph trajectory
drawImage(handles);
drawGraph(handles)


function handles = loadImage(handles)

% 1) Get time/XY (i/j) indicies for new image + label matricies
i = handles.xy;
xy_idx = find(handles.xys==handles.xy,1,'first');
j = handles.time-handles.shift(xy_idx);
img_path = [handles.locations.scope, filesep, handles.parameters.ImagePath,filesep,...
    eval(handles.ImageExpr)];
% 2) Get bit depth of image, load, and rotate (if necessary)
if ~isfield(handles.parameters,'BitDepth')
    imfo = imfinfo(img_path);
    handles.parameters.BitDepth = imfo.BitDepth;
end
measure_img = checkread(img_path,handles.parameters.BitDepth);
if size(measure_img,1)>size(measure_img,2)
    measure_img = imrotate(measure_img,90);
    flip_flag = 1;
end
% 3) Load label matricies (rotate if necessary)
load([handles.OutputDir,filesep,'xy',num2str(i),filesep,...
    'NuclearLabels',filesep,'NuclearLabel-',numseq(j,4),'.mat'])
try
    load([handles.OutputDir,'xy',num2str(i),filesep,...
        'CellLabels',filesep,'CellLabel-',numseq(j,4),'.mat'])
catch me
    CellLabel = NuclearLabel;
end
if exist('flip_flag','var') && (flip_flag ==1)
    handles.nuc_label = imrotate(NuclearLabel,90);
    handles.cell_label = imrotate(CellLabel,90);
else
    handles.nuc_label = NuclearLabel;
    handles.cell_label = CellLabel;
end

% 4) Alter base image for display: saturate according to an early + a late image
if ~isfield(handles,'imgmax')
    j = handles.parameters.TimeRange(round(end/2));
    img_path = namecheck([handles.locations.scope,'/', handles.parameters.ImagePath,'/',eval(handles.ImageExpr)]);
    scale_img = checkread(img_path,handles.parameters.BitDepth);
    rng1 = prctile(scale_img(:),[3 99]);
    j = handles.parameters.TimeRange(1);
    img_path = namecheck([handles.locations.scope,'/', handles.parameters.ImagePath,'/',eval(handles.ImageExpr)]);
    scale_img = checkread(img_path,handles.parameters.BitDepth);
    rng2 = prctile(scale_img(:),[3 99]);
    handles.imgmax = max([rng1 rng2]);
    handles.imgmin = min([rng1 rng2]);  
end
measure_img = (measure_img-handles.imgmin)/(handles.imgmax-handles.imgmin);
measure_img(measure_img>1) = 1; measure_img(measure_img<0) = 0;
handles.img = uint8(measure_img*255);

% 5) Get centroids and labels
tmp_label = double(handles.nuc_label);
tmp_label(~ismember(tmp_label,handles.cells_all)) = 0;
props = regionprops(tmp_label,'Centroid');
tmpcell = struct2cell(props);
tmpmat = cell2mat(tmpcell(1,:));
handles.centroids = [tmpmat(1:2:end)', tmpmat(2:2:end)'];
handles.text_labels = cellstr(num2str(unique(tmp_label(tmp_label>0))));

function drawImage(handles)
% 1) Make overlays for cell/nuclear labels
%%
cell_all = ismember(handles.cell_label,handles.cells_all);
nuc_all = ismember(handles.nuc_label,handles.cells_all);

cell_selected = ismember(handles.cell_label,unique(handles.lineage(handles.xy_subset(handles.cell_idx))));
nuc_selected = ismember(handles.nuc_label,unique(handles.lineage(handles.xy_subset(handles.cell_idx))));
%%
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

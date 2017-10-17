function varargout = checkDynamics(varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = checkDynamics(varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKDYNAMICS creates a figure showing tracked/outlined cells that corrrespond to a 
% dynamic measurement, graphed below. Measurements are first processed by specific 
% vizualization functions called (e.g.) see_nfkb or see_fret.
%
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

% 3) Load specified measurements file. Get source image directory and output label matricies
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


% 4) Load visualization data into handles, set GUI elements
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

if isfield(handles,'imgmax')
    handles = rmfield(handles,'imgmax');
end
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
handles.see_fcn = fcn_name;


% Get primary fields- timepoints & data
handles.t_img = info.parameters.TimeRange;
handles.t_graph = graph.t;

% Handle multi-timepoint measurement (if applicable)
if isfield(info,'t_aux'); 
    handles.t_aux = info.t_aux;
else
    handles.t_aux = handles.t_img;
end

handles.xys = unique(graph.celldata(:,1));
handles.celldata = graph.celldata;
handles.var = graph.var;


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
    handles.ylim = prctile(handles.var(:),[1 99]);
end

% Make calculated fields - average across all cells per timept
handles.mu = nanmean(handles.var);
handles.sigma = nanstd(handles.var);
handles.graph_lim = prctile(handles.var,[1 99.5]);

c = loadcolormaps;
handles.cmap = c.red_hot;

% Initialize slider + popup values
if length(handles.xys)>1
    set(handles.slider1, 'Min',1, 'Max',length(handles.xys),'SliderStep',[1 4]/length(handles.xys),'Value',1);
end
handles.xy = handles.xys(1);
set(handles.text5,'String',['xy ' ,num2str(handles.xy)]);
set(handles.text7,'String',['t = ' ,num2str(min(handles.t_img))]);

v1 = 1; v2 = length(handles.t_graph); 
set(handles.slider2, 'Min',v1, 'Max',v2,'SliderStep',[1/(v2-v1) 4/(v2-v1)],'Value',v1);
handles.timept = v1;
% Initialize image+graph, update handles
handles = newXY(handles);
handles_out = handles;
guidata(handles.figure1, handles)



% - - - - - - GUI functions: update values, call appropriate drawing functions - - - - - - - 
function slider1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderval = round(get(hObject,'Value'));
handles.xy = handles.xys(sliderval);
set(hObject,'Value',sliderval);
set(handles.text5,'String',['xy ' , num2str(handles.xy)]);
handles = newXY(handles);
guidata(handles.figure1, handles);


function slider2_Callback(hObject, eventdata, handles)
handles.timept = round(get(hObject,'Value'));
set(handles.text7,'String',['t = ',num2str(handles.t_aux(handles.timept))]);
set(hObject,'Value',handles.timept);
handles = newTime(handles);
guidata(handles.figure1, handles);

function popupmenu1_Callback(hObject, eventdata, handles)
cellno = get(handles.popupmenu1,'Value');
if cellno ~= handles.cell_idx
    handles.cell_idx = cellno;
    handles.cell = unique(handles.lineage(handles.xy_subset(handles.cell_idx),:));
    handles.cell(isnan(handles.cell)) = [];
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


% - - - - - IMAGE/GRAPH UPDATE functions: - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
handles.cell_all = unique(handles.lineage(handles.xy_subset,:));
handles.cell_all(isnan(handles.cell_all)) = [];
handles.cell = unique(handles.lineage(handles.xy_subset(1),:));
handles.cell(isnan(handles.cell)) = [];
set(handles.popupmenu1,'Value',1);

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
xy_idx = find(handles.xys==handles.xy,1,'first'); % XY index
j= handles.timept-handles.shift(xy_idx); % Timepoint index (for images)
label_idx = find( abs(handles.t_img-j) == min(abs(handles.t_img-j))); % Timepoint index (for label mats)

% For multi-TP need to distinguish btw the timepoint for image vs timept for mask.
% For now, we're going to assume shift isn't compatible w/ this - would require some jankiness otherwise.
if ~isequal(handles.t_img, handles.t_aux)
    j_aux = handles.t_aux(handles.timept);
    j = floor(j_aux);
else
    j_aux = j;
end

% Load corresponding image (defined by i/j pair). If the viz function ("see_...") returns a basic ImageExpr, load this.
% If it doesn't (e.g. because the image must be derived, like a FRET image), run the viz function in 'GetImage' mode.

if ~isempty(handles.ImageExpr)
    img_path = namecheck([handles.locations.scope, filesep, handles.parameters.ImagePath,filesep,...
        eval(handles.ImageExpr)]);
    if ~isfield(handles.parameters,'BitDepth')
        imfo = imfinfo(img_path);
        handles.parameters.BitDepth = imfo.BitDepth;
    end
    measure_img = checkread(img_path,handles.parameters.BitDepth);
    
    if ~isfield(handles,'imgmax')
        j = handles.t_img(round(end/2));
        img_path = namecheck([handles.locations.scope,'/', handles.parameters.ImagePath,'/',eval(handles.ImageExpr)]);
        scale_img = checkread(img_path,handles.parameters.BitDepth);
        rng1 = prctile(scale_img(:),[3 99]);
        j = handles.t_img(1);
        img_path = namecheck([handles.locations.scope,'/', handles.parameters.ImagePath,'/',eval(handles.ImageExpr)]);
        scale_img = checkread(img_path,handles.parameters.BitDepth);
        rng2 = prctile(scale_img(:),[3 99]);
        handles.imgmax = max([rng1 rng2]);
        handles.imgmin = min([rng1 rng2]);  
    end
  
else 
    measure_img = eval([handles.see_fcn,'(handles.AllMeasurements, ''GetImage'', [i j_aux]);']);  
    if ~isfield(handles,'imgmax')
        j_aux = handles.t_img(round(end/2));
        scale_img = eval([handles.see_fcn,'(handles.AllMeasurements, ''GetImage'', [i j_aux]);']);
        rng1 = prctile(scale_img(:),[3 99]);
        j_aux = handles.t_img(1);
        scale_img = eval([handles.see_fcn,'(handles.AllMeasurements, ''GetImage'', [i j_aux]);']);
        rng2 = prctile(scale_img(:),[3 99]);
        handles.imgmax = max([rng1 rng2]);
        handles.imgmin = min([rng1 rng2]);
    end
end


if size(measure_img,1)>size(measure_img,2)
    measure_img = imrotate(measure_img,90);
    flip_flag = 1;
end
% Resize to a max size of ~1024 px width (our max plotting size)
if size(measure_img,2)>1200
    scale_factor = 1024/size(measure_img,2);
    measure_img = imresize(measure_img,scale_factor);
else
    scale_factor = 1;
end

% 3) Load label matricies (rotate if necessary)
load([handles.OutputDir,filesep,'xy',num2str(i),filesep,...
    'NuclearLabels',filesep,'NuclearLabel-',numseq(label_idx,4),'.mat'])
try
    load([handles.OutputDir,'xy',num2str(i),filesep,...
        'CellLabels',filesep,'CellLabel-',numseq(label_idx,4),'.mat'])
catch me
    CellLabel = [];
end

% 4) Process & scale image for display
if exist('flip_flag','var') && (flip_flag ==1)
    handles.nuc_label = imrotate(NuclearLabel,90);
    handles.cell_label = imrotate(CellLabel,90);
else
    handles.nuc_label = NuclearLabel;
    handles.cell_label = CellLabel;
end

if scale_factor <1
    handles.nuc_label = imresize(handles.nuc_label,scale_factor,'nearest');
    if ~isempty(handles.cell_label)
        handles.cell_label = imresize(handles.cell_label,scale_factor,'nearest');
    end
end
measure_img = (measure_img-handles.imgmin)/(handles.imgmax-handles.imgmin);
measure_img(measure_img>1) = 1; measure_img(measure_img<0) = 0;
handles.img = uint8(measure_img*255);

% 5) Get centroids and labels
tmp_label = double(handles.nuc_label);
tmp_label(~ismember(tmp_label,handles.cell_all)) = 0;
props = regionprops(tmp_label,'Centroid');
tmpcell = struct2cell(props);
tmpmat = cell2mat(tmpcell(1,:));
handles.centroids = [tmpmat(1:2:end)', tmpmat(2:2:end)'];
handles.centroids(isnan(handles.centroids(:,1)),:) = [];
handles.text_labels = cellstr(num2str(unique(tmp_label(tmp_label>0))));

function drawImage(handles)
% 1) Create masks for active cells/nuclei
if ~isempty(handles.cell_label)
    tmp1 = handles.cell_label;
    tmp1(~ismember(tmp1,handles.cell_all)) = 0;
    tmp1(ismember(tmp1,handles.cell)) = 0;
    borders_cell = (imdilate(tmp1,ones(3))-tmp1)>0;
end

tmp1 = handles.nuc_label;
tmp1(~ismember(tmp1,handles.cell_all)) = 0;
tmp1(ismember(tmp1,handles.cell)) = 0;
borders_nuc = (imdilate(tmp1,ones(3))-tmp1)>0;
nucmask =  ismember(handles.nuc_label,handles.cell);
if isempty(handles.cell_label)
    borders_individ = (imdilate(nucmask,ones(3))&~nucmask);  
else
    cellmask = ismember(handles.cell_label,handles.cell);
    borders_individ = (imdilate(cellmask,ones(3))&~cellmask) | (imdilate(nucmask,ones(3))&~nucmask);    
end

if ~isempty(get(handles.axes1,'Children'))
    xlim = get(handles.axes1,'XLim');
    ylim = get(handles.axes1,'YLim');
    setflag = 1;
else
    setflag = 0;  
end

% 2) Overlay borders + mask of selected cell
% [ NEON theme: blue: [0 230 240] | yellow - [255 255 0] | green - [0 255 0] | orange - [255 104 0] ]
neon1 = [0 230 240]; 
neon2 = [255 255 0];
neon3 = [255 255 255];

RGB = maskoverlay(handles.img, imdilate(borders_nuc,ones(3)), neon2,0.15);
RGB = maskoverlay(RGB, borders_nuc, neon2, 0.4);
% Overlay nuclear labels (idf distinct from cells)
if ~isempty(handles.cell_label)

    % Overlay borders for (filtered) cell boundaries
    RGB = maskoverlay(RGB, imdilate(borders_cell,ones(3)), neon1, 0.15);
    RGB = maskoverlay(RGB, borders_cell, neon1, 0.4);
end
% Make object region colormapped
RGB2 = uint8(cat(3,reshape(handles.cmap(handles.img+1,1),size(handles.img)),...
    reshape(handles.cmap(handles.img+1,2),size(handles.img)),...
    reshape(handles.cmap(handles.img+1,3),size(handles.img)))*255);
RGB = RGB.*(repmat(uint8(~cellmask),[1 1 3])) +  RGB2.*(repmat(uint8(cellmask),[1 1 3]));
RGB = maskoverlay(RGB,borders_individ,neon3, 0.6);

% 3) Add text, draw image
RGB = insertText(RGB, handles.centroids,handles.text_labels,'TextColor',[247 197 113],'BoxOpacity',0,...
    'AnchorPoint','Center','Font','Arial','FontSize',16);
imshow(RGB,'Parent',handles.axes1)
set(handles.axes1,'XTick',[],'YTick',[])
if setflag
    set(handles.axes1,'xlim',xlim,'ylim',ylim,'XTick',[],'YTick',[])
else
    axis image
end

% text(handles.centroids(:,1),handles.centroids(:,2),handles.text_labels,'Color',[1 1 1],...
%     'FontName','Arial Narrow','FontSize',14,'Parent',handles.axes1);

function drawGraph(handles)
% 1) Plot field-of-view subpopulation as area plot
clr = [0.6588 0.7059 0.8000];
fill_clr = (clr+3*[0.9 0.9 0.9])/4;

cla(handles.axes2)
plot(handles.t_graph,handles.mu_subset,'Color',clr,'LineWidth',2,'Parent',handles.axes2)
hold(handles.axes2,'on')

patch([handles.t_graph,fliplr(handles.t_graph)],...
    [handles.mu_subset+handles.sigma_subset,fliplr(handles.mu_subset-handles.sigma_subset)],...
    fill_clr,'LineStyle','none','Parent',handles.axes2)
plot(handles.t_graph,handles.mu_subset,'Color',clr,'LineWidth',2,'Parent',handles.axes2)

% 2) Plot whole-population average as dashed/dotted lines
plot(handles.t_graph,handles.mu,'Color',[.2 .2 .2],'LineStyle','-.','Parent',handles.axes2)
% plot(handles.t_graph,handles.mu+handles.sigma,'Color',[.2 .2 .2],'LineStyle','-.','Parent',handles.axes2)
% plot(handles.t_graph,handles.mu-handles.sigma,'Color',[.2 .2 .2],'LineStyle','-.','Parent',handles.axes2)

% 3) Plot single cell as solid line 
plot(handles.axes2,handles.t_graph, handles.var(handles.xy_subset(handles.cell_idx),:),'Color',[68 129 227]/255,'LineWidth',3)


% Make a line for timept
plot(handles.t_graph(handles.timept),handles.var(handles.xy_subset(handles.cell_idx),handles.timept),'o',...
    'MarkerSize',16,'Color','k','Parent',handles.axes2)
handles.line1 = line([handles.t_graph(handles.timept) handles.t_graph(handles.timept)],handles.ylim,...
    'Color',[0 0 0],'Parent',handles.axes2);
hold(handles.axes2,'off')

set(handles.axes2,'FontSize',12, 'YLim',handles.ylim,'XLim',[min(handles.t_graph),max(handles.t_graph)],'Box','on')

text(max(handles.t_graph),max(handles.ylim),...
    ['x =  ',num2str(.01*round(handles.var(handles.xy_subset(handles.cell_idx),handles.timept)*100)),' '],...
    'VerticalAlignment','top','HorizontalAlignment','right','Parent',handles.axes2,'FontSize',16,'FontWeight','bold',...
    'Color',[68 129 227]/255)

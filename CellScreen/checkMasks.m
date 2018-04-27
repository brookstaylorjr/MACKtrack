function [] = checkMasks(track_folder, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = showmasks(segmentation_folder)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% CHECKMASKS is akin to checkDynamics, but is specific to "screen" functionality of MACKtrack.
% It plots nuclear/cytoplasmic boundaries on original images, enabling user to check their quality. 
% 
% INPUTS
% track_folder folder containing AllData structure, masks, etc
% parameters   tracking parameters - can be found (as 'TrackingParameters.mat') in same directory as AllData
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Load in parameters, if not provided
if nargin<2
    disp(['Loading in "TrackingParameters.mat" from ', track_folder])
    load(namecheck([track_folder,filesep,'TrackingParameters.mat']))
end


% Load scope/data locations
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end
handles.data_dir = track_folder;

%% FOR TESTING ONLY
% track_folder = '/Volumes/labdata/brooks/Tracked/2016-10-24_physio-screen/Day 1/';
% tmp = load(namecheck([track_folder,filesep,'TrackingParameters.mat']));
% parameters = tmp.parameters;
% load locations.mat
% parameters.ImagePath_full = parameters.ImagePath_full(13:end); % tmp bugfix -> get frid of /mnt/fluffy.



%% Get list of all sites (from 'NuclearLabels') - assume MicroXL formatting (e.g. '_B01_s2')
handles.labels_n = quickdir(namecheck([track_folder,filesep,'NuclearLabels']));
handles.labels_n(cellfun(@isempty,regexp(handles.labels_n,'NuclearLabel.*.mat'))) = [];
if exist([track_folder,filesep,'CellLabels'],'dir')
   handles.labels_c = quickdir(namecheck([track_folder,filesep,'CellLabels']));
   handles.labels_c(cellfun(@isempty,regexp(handles.labels_c,'CellLabel.*.mat'))) = []; 
end
%%
get_substr = @(str1) str1(14:end-4);
handles.wells = cellfun(get_substr, handles.labels_n,'UniformOutput',0);


handles.image_dir = namecheck([locations.scope, filesep, parameters.ImagePath_full]);



image_list = quickdir(handles.image_dir);
image_list(cellfun(@isempty,strfind(image_list,'.tif'))) = []; % Drop non-tif images

[image_list,scope] = wellmatch(image_list,'');
%%
if strcmpi(scope,'metaxpress')
    conv_well = @(str1) ['_',str1(1:3),'_s',num2str(eval(str1(end-1:end)))];
elseif strcmpi(scope,'slidebook')
    conv_well = @(str1) ['- ',str1(1),num2str(eval(str1(2:3))),'_',numseq(-1+eval(str1(end-1:end)),3)];
else
    error('Unsupported naming type found for your images. See ''wellmatch'' function for more details')

end
handles.wells = cellfun(conv_well,handles.wells,'UniformOutput',0);



%%
% Get list of images that correspond to each state (omit thumbnails)
handles.images = cell(size(handles.wells));
for i = 1:length(handles.wells)
    handles.images{i} = image_list(~cellfun(@isempty,strfind(image_list,handles.wells{i})) & ...
        cellfun(@isempty,strfind(image_list,'_thumb')));
end

%%

colormaps = loadcolormaps;
handles.cmap = colormaps.purple_hot;
handles.radius = parameters.MinNucleusRadius;

%%  Create figure, axes, and GUI elements
zlength = length(handles.wells);
handles.figure1 = figure('Position',positionfig(800,600));
handles.axes1 = axes('Parent',handles.figure1);

% Create status text
handles.text1 = uicontrol('Style','text','String','','BackgroundColor',[1 1 1]);

% Create GUI elements
handles.slider1  = uicontrol('Style', 'slider','Max',zlength,'Min',1,'Value',1,'SliderStep',[1/(zlength-1) 1/(zlength-1)],'BackgroundColor',[.99 .99 1]);
set(handles.text1,'String',[handles.wells{1}])
handles.popup1  = uicontrol('Style', 'popup','Value',1,'String',handles.images{1},'BackgroundColor',[.99 .99 1]);
handles.popup2  = uicontrol('Style', 'popup','Value',1,'String',{'Nucleus', 'Cell', 'Annulus'} ,'BackgroundColor',[.99 .99 1]);
handles.reset  = uicontrol('Style', 'pushbutton','String','Reset zoom','BackgroundColor',[.99 .99 1]);

% Create image
[handles.xlim, handles.ylim] = displayImage(handles);


% Create all callbacks
set(handles.popup1,'Callback',{@popup1_Callback,handles});
set(handles.popup2,'Callback',{@popup2_Callback,handles});
set(handles.slider1,'Callback',{@slider_Callback,handles});
set(handles.reset,'Callback',{@reset_Callback,handles});
set(handles.figure1,'ResizeFcn',{@fig_resize,handles},'Toolbar','figure');


% Trigger resize element
fig_resize([],[],handles)




function slider_Callback(hObj,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update contents of figure based on slider's position
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
wellVal = round(get(hObj,'Value'));
set(hObj,'Value',wellVal)

set(handles.text1,'String',[handles.wells{wellVal}])
set(handles.popup1,'String',handles.images{wellVal});


xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
displayImage(handles, xlim, ylim);
% ========================================================================================



function popup1_Callback(hObj,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update contents of figure based on slider's position
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
displayImage(handles, xlim, ylim);

% ========================================================================================


function popup2_Callback(hObj,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update contents of figure based on slider's position
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
displayImage(handles, xlim, ylim);

% ========================================================================================

function reset_Callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Reset zoom state of figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xlim = handles.xlim;
ylim = handles.ylim;
set(gca,'xlim',xlim,'ylim',ylim)

% ========================================================================================


function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize figure and  slider smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.figure1,'Position');
h = figPos(4);
w = figPos(3);
w_ref = round(w/2) - 410;

set(handles.slider1,'Position',[w_ref 2 200 20])
set(handles.text1,'Position',[w_ref+210 2 60 20]);
set(handles.popup1,'Position',[w_ref+280 5 300 20])
set(handles.popup2,'Position',[w_ref+590 5 100 20])
set(handles.reset,'Position',[w_ref+710 8 66 18]);


axis image
set(handles.axes1,'OuterPosition',[0 20/h 1 (h-32)/h],'LooseInset',get(handles.axes1,'TightInset')+[10/w 10/h 10/w 10/h])
set(handles.axes1,'YTick',[],'XTick',[])
% ========================================================================================


function  [xlim, ylim] = displayImage(handles, xlim, ylim)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Display image (with specified boundaries overlaid)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
wellVal = get(handles.slider1,'Value');
imgVal = get(handles.popup1,'Value');
boundaryVal = get(handles.popup2,'Value');

% Load constituent parts of image - image 'base', Nuclear and Cell labels (if present)
base_img = checkread(namecheck([handles.image_dir,filesep,handles.images{wellVal}{imgVal}]));
if nargin<3
    xlim = [1 size(base_img,2)];
    ylim = [1 size(base_img,1)];
end


if ~isfield(handles,'labels_c')
    boundaryVal(boundaryVal==2) = 3;
end

switch boundaryVal
    case 1 % Nuclear boundaries
        load(namecheck([handles.data_dir,filesep,'NuclearLabels',filesep,handles.labels_n{wellVal}]))
        boundaries = NuclearLabel;
    case 2 % Cell boundaries
        load(namecheck([handles.data_dir,filesep,'CellLabels',filesep,handles.labels_c{wellVal}]))
        boundaries = CellLabel;
    case 3 % Annulus boundaries
        load(namecheck([handles.data_dir,filesep,'NuclearLabels',filesep,handles.labels_n{wellVal}]))
        % Create annulus to estimate cytoplasmic area.
        tmp_mask = imdilate(NuclearLabel>0,diskstrel(handles.radius*1.5)); % (Expand by 1.5 radii)
        boundaries = IdentifySecPropagateSubfunction(double(NuclearLabel),zeros(size(NuclearLabel)),tmp_mask,100);
end
boundaries = (imdilate(boundaries,ones(3))-boundaries)>0;


% Map orig image to a colormap (purple_hot)
bounds = prctile(base_img(:),[3 98]);

base_img = (base_img-min(bounds))/range(bounds);
base_img(base_img>11) = 1; base_img(base_img<0) = 0;
base_img = uint8(base_img*255);

disp_image = uint8(cat(3,reshape(handles.cmap(base_img+1,1),size(base_img)),...
    reshape(handles.cmap(base_img+1,2),size(base_img)),...
    reshape(handles.cmap(base_img+1,3),size(base_img)))*255);
disp_image = maskoverlay(disp_image, boundaries, [248   152    29], 0.6);

imshow(disp_image,'Parent',handles.axes1)
axis image
set(gca,'xlim',xlim,'ylim',ylim)
set(handles.axes1,'YTick',[],'XTick',[])

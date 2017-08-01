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
    disp('Loading in default "TrackingParameters.mat"')
    load(namecheck([track_folder,filesep,'TrackingParameters.mat']))
end







%% 

% Get parameters for this tracking set
tmp = load(namecheck([track_folder,filesep,'TrackingParameters.mat']));
parameters = tmp.parameters;

% Load scope/data locations
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end


%% FOR TESTING ONLY
track_folder = '/Volumes/labdata/brooks/Tracked/2016-10-24_physio-screen/Day 1/';
tmp = load(namecheck([track_folder,filesep,'TrackingParameters.mat']));
parameters = tmp.parameters;
load locations.mat
parameters.ImagePath_full = parameters.ImagePath_full(13:end); % tmp bugfix -> get frid of /mnt/fluffy.



%% Get list of all sites (from 'NuclearLabels') - assume MicroXL formatting (e.g. '_B01_s2')
handles.labels_n = quickdir(namecheck([track_folder,filesep,'NuclearLabels']));
handles.labels_n(cellfun(@isempty,regexp(handles.labels_n,'NuclearLabel.*.mat'))) = [];
if exist([track_folder,filesep,'CellLabels'],'dir')
   handles.labels_c = quickdir(namecheck([track_folder,filesep,'CellLabels']));
   handles.labels_c(cellfun(@isempty,regexp(handles.labels_c,'CellLabel.*.mat'))) = []; 
end

get_substr = @(str1) str1(14:end-4);
handles.wells = cellfun(get_substr, label_list,'UniformOutput',0);
conv_well = @(str1) ['_',str1(1:3),'_s',num2str(eval(str1(end-1:end)))];
handles.wells = cellfun(conv_well,handles.wells,'UniformOutput',0);

% Get list of images that correspond to each state (omit thumbnails)
handles.image_dir = namecheck([locations.scope, filesep, parameters.ImagePath_full]);
image_list = quickdir(handles.image_dir);
image_list(cellfun(@isempty,strfind(image_list,'.tif'))) = []; % Drop non-tif images
handles.images = cell(size(handles.wells));
for i = 1:length(handles.wells)
    handles.images{i} = image_list(~cellfun(@isempty,strfind(image_list,handles.wells{i})) & ...
        cellfun(@isempty,strfind(image_list,'_thumb')));
end

%%
handles.wellVal = 1;
handles.imgVal = 1;
handles.boundaryVal = 1;
colormaps = loadcolormaps;
handles.cmap = colormaps.purple_hot;







%%  Create figure, axes, and GUI elements
zlength = length(handles.wells);

% Load 1st image of 1st well
handles.image_curr = loadImage(handles);

%%



handles.figure1 = gcf;
handles.axes1 = axes('Parent',handles.figure1);
% Create status text
handles.text1 = uicontrol('Style','text','String','','BackgroundColor',[1 1 1]);
% Create slider
handles.slider1  = uicontrol('Style', 'slider','Max',zlength,'Min',1,'Value',1,'SliderStep',[1/(zlength-1) 1/(zlength-1)],'BackgroundColor',[.99 .99 1]);
set(handles.slider1,'Callback',{@slider_Callback,handles});
% Make "Reset zoom" button
handles.reset  = uicontrol('Style', 'pushbutton','String','Reset zoom','BackgroundColor',[.99 .99 1]);
set(handles.reset,'Callback',{@reset_Callback,handles});





function slider_Callback(hSlide,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update contents of figure based on slider's position
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.wellVal = round(get(hSlide,'Value'));
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
displayImage(handles)
set(handles.text1,'String',[handles.wells{sliderVal},'/',num2str(size(handles.input,3))])
axis image
set(gca,'xlim',xlim,'ylim',ylim)
set(handles.axes1,'YTick',[],'XTick',[])
% ========================================================================================

function dropdown1_Callback(hSlide,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update contents of figure based on slider's position
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles.wellVal = round(get(hSlide,'Value'));
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
displayImage(handles)
set(handles.text1,'String',[handles.wells{sliderVal},'/',num2str(size(handles.input,3))])
axis image
set(gca,'xlim',xlim,'ylim',ylim)
set(handles.axes1,'YTick',[],'XTick',[])
% ========================================================================================











function reset_Callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Reset zoom state of figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

xlim = [1 size(handles.image_curr,2)];
ylim = [1 size(handles.image_curr,1)];

set(gca,'xlim',xlim,'ylim',ylim)

% ========================================================================================


function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize figure and  slider smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.figure1,'Position');
h = figPos(4);
w = figPos(3);
set(handles.slider1,'Position',[floor(w/2)-300 5 400 20])
set(handles.reset,'Position',[floor(w/2)+182,8, 60, 18]);
set(handles.text1,'Position',[floor(w/2)+102,5, 80, 20]);
axis image
set(handles.axes1,'OuterPosition',[0 20/h 1 (h-32)/h],'LooseInset',get(handles.axes1,'TightInset')+[10/w 10/h 10/w 10/h])
set(handles.axes1,'YTick',[],'XTick',[])
% ========================================================================================


function  disp_image = displayImage(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Display image (with specified boundaries overlaid)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Load constituent parts of image - image 'base', Nuclear and Cell labels (if present)
base_img = checkread(namecheck([handles.image_dir,filesep,handles.images{handles.wellVal}{handles.imgVal}]));

boundaryVal = handles.boundaryVal;
if ~isfield(handles,'labels_c')
    boundaryVal(boundaryVal==2) = 3;
end

switch boundaryVal
    case 1 % Nuclear boundaries
        load(namecheck([track_folder,filesep,'NuclearLabels',filesep,handles.labels_n{handles.wellVal}]))
        boundaries = NuclearLabel;
    case 2 % Cell boundaries
        load(namecheck([track_folder,filesep,'CellLabels',filesep,handles.labels_c{handles.wellVal}]))
    case 3 % Annulus boundaries
        load(namecheck([track_folder,filesep,'NuclearLabels',filesep,handles.labels_n{handles.wellVal}]))
        % Create annulus to estimate cytoplasmic area.
        tmp_mask = imdilate(NuclearLabel>0,diskstrel(parameters.MinNucleusRadius*1.5)); % (Expand by 1.5 radii)
        boundaries = IdentifySecPropagateSubfunction(double(lNuclearLabel),zeros(size(lNuclearLabel)),tmp_mask,100);
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


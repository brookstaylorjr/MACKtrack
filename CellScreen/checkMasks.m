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




%% Get list of all sites (from 'NuclearLabels') - assume MicroXL formatting (e.g. '_B01_s2')
label_list = quickdir(namecheck([track_folder,filesep,'NuclearLabels']));
label_list(cellfun(@isempty,regexp(label_list,'NuclearLabel.*.mat'))) = [];
if exist([track_folder,filesep,'CellLabels'],'dir')
   label_list = quickdir(namecheck([track_folder,filesep,'CellLabels']));
   label_list(cellfun(@isempty,regexp(label_list,'CellLabel.*.mat'))) = []; 
end

get_substr = @(str1) str1(14:end-4);
handles.wells = cellfun(get_substr, label_list,'UniformOutput',0);
conv_well = @(str1) ['_',str1(1:3),'_s',num2str(eval(str1(end-1:end)))];
handles.wells = cellfun(conv_well,handles.wells,'UniformOutput',0);

%% Get list of images that correspond to each state.
handles.image_dir = namecheck([locations.scope, filesep, parameters.ImagePath_full]);
image_list = quickdir(handles.image_dir);
image_list(cellfun(@isempty,strfind(image_list,'.tif'))) = []; % Drop non-tif images
handles.images = cell(size(handles.wells));
for i = 1:length(handles.wells)
    handles.images{i} = image_list(~cellfun(@isempty,strfind(image_list,handles.wells{i})));
end


%%  Create figure, axes, and GUI elements
zlength = length(handles.wells);

% Load 1st image of 1st well
handles.image_curr = checkread(



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
sliderVal = round(get(hSlide,'Value'));
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
if handles.flag1
    imagesc(handles.input(:,:,sliderVal),'Parent',handles.axes1,handles.bounds)
else
    imagesc(handles.input(:,:,sliderVal),'Parent',handles.axes1)
end
set(handles.text1,'String',[handles.wells{sliderVal},'/',num2str(size(handles.input,3))])
axis image
set(gca,'xlim',xlim,'ylim',ylim)
set(handles.axes1,'YTick',[],'XTick',[])


% ========================================================================================

function reset_Callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Reset zoom state of figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

xlim = [1 size(handles.input,2)];
ylim = [1 size(handles.input,1)];

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

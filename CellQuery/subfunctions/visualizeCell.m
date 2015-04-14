function  visualizeCell(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% visualizeCell(handles)
% VISUALIZECELL makes new figure to show/cycle outlined cell on user's choice of images
%
% INPUT:
% handles      main data structure provided by UCSDcellQuery GUI
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Pull off necessary data from handles, and create a new figure
visualize.Options = handles.Options;
visualize.parameters = handles.trackparameters;
visualize.image_dir = [handles.locations.scope,visualize.parameters.ImagePath];

visualize.figure1 = figure;
set(visualize.figure1,'Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')

% Find directory where masks are saved
if ~ismember(visualize.Options.VisualizeXY,unique(handles.Measurements.CellData(:,1)))
    error('The specified XY was not found in this experiment.')
end
slash_idx = strfind(handles.file_loc,filesep);
if ~exist([handles.file_loc(1:slash_idx(end)),filesep,'xy',num2str(visualize.Options.VisualizeXY)],'dir')
    visualize.label_dir = [handles.locations.data,visualize.parameters.SaveDirectory];
else
    visualize.label_dir = handles.file_loc(1:slash_idx(end));
end
if ~exist(visualize.label_dir,'dir')
    error(['Unable to find "', visualize.label_dir,'"'])
end

% Initialize i/j to construct image names
i = min(handles.trackparameters.TimeRange);
j = handles.Options.VisualizeXY;

% Create the figure axes
visualize.axes1 = axes('Parent',visualize.figure1);


% Create slider and set step size to correspond to number of frames andalizyed
num_frames = length(handles.trackparameters.TimeRange) - handles.trackparameters.StackSize+1;
visualize.slider1  = uicontrol('Style', 'slider','Max',num_frames,'Min',1,'Value',1,'SliderStep',[1/(num_frames-1) 10/(num_frames-1)],'BackgroundColor',[.99 .99 1]);

% Create dropdown and initialize with Nucleus, Phase, and any defined Aux1 images
visualize.imageList = {'Nucleus', 'Cell'};
modules = visualize.parameters.ModuleNames;
for ind1 = 1:length(modules)
    % If the module is in use, and auxImageExpr is a file add it to the list of images
    if handles.trackparameters.(modules{ind1}).Use 
        try
            if exist([visualize.image_dir,filesep,eval(visualize.parameters.(modules{ind1}).ImageExpr)],'file')
                visualize.imageList = cat(2,visualize.imageList,[modules{ind1},'Img']);
            end
        catch me
        end
    end
end

visualize.popup1  = uicontrol('Style', 'popupmenu','String',visualize.imageList,'FontName','sans serif',...
    'BackgroundColor',[.99 .99 1]);

% Set properties and callbacks, then run initial display
set(visualize.axes1,'YTick',[],'XTick',[])
set(visualize.slider1,'Callback',{@ui_Callback,visualize});
set(visualize.popup1,'Callback',{@ui_Callback,visualize});
set(visualize.figure1,'ResizeFcn',{@fig_resize,visualize},'Toolbar','figure');

% Before first call, grab cell's CellData and create axes title
notfound = 0;
xyset = unique(handles.Measurements.CellData(:,1));
if ~ismember(visualize.Options.VisualizeXY,xyset)
    visualize.Options.VisualizeXY = min(xyset);
    notfound = 1;
end
row1 = find((handles.Measurements.CellData(:,1)==visualize.Options.VisualizeXY) & ...
    (handles.Measurements.CellData(:,2)==visualize.Options.VisualizeCell),1,'first');
if isempty(row1)
    visualize.Options.VisualizeCell = ...
        min(handles.Measurements.CellData(handles.Measurements.CellData(:,1)==visualize.Options.VisualizeXY,2));
    row1 = find((handles.Measurements.CellData(:,1)==visualize.Options.VisualizeXY) &...
        (handles.Measurements.CellData(:,2)==visualize.Options.VisualizeCell),1,'first');
    notfound = 1;
end
visualize.title = ['frame #',numseq(1,3),', xy #', num2str(visualize.Options.VisualizeXY),...
    ', cell #', num2str(visualize.Options.VisualizeCell),...
    ' [in frames ',num2str(handles.Measurements.CellData(row1,3)),'-',num2str(handles.Measurements.CellData(row1,4)),']'];
if notfound
    visualize.title = [visualize.title, '(selected cell/XY not found)'];
end

draw_figure(visualize,1,visualize.imageList{1})
fig_resize([],[],visualize)


function ui_Callback(~,~,visualize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Read slider position and popup string, update figure accordingly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
frameNo = round(get(visualize.slider1,'Value'));
imageStr = visualize.imageList{get(visualize.popup1,'Value')};
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
draw_figure(visualize, frameNo, imageStr)
set(gca,'xlim',xlim,'ylim',ylim)

% ========================================================================================


function fig_resize(~,~,visualize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize figure,  slider  and dropdown menu smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(visualize.figure1,'Position');
h = figPos(4);
w = figPos(3);
set(visualize.slider1,'Position',[floor(w/2)-200 5 200 20])
set(visualize.popup1,'Position',[floor(w/2)+10 5 190 20])
axis image
set(visualize.axes1,'OuterPosition',[0 20/h 1 (h-32)/h],'LooseInset',get(visualize.axes1,'TightInset')+[10/w 10/h 10/w 10/h])
set(visualize.axes1,'YTick',[],'XTick',[])
% ========================================================================================


function draw_figure(visualize, frameNo, imageStr)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Display image based on frameNo and 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
i = visualize.Options.VisualizeXY;
j = frameNo;
cellno = visualize.Options.VisualizeCell;

% Grab old title and update frame number
if ~isfield(visualize,'title')
    visualize.title = get(get(visualize.axes1,'Title'),'String');
end
visualize.title(8:10) = numseq(j,3);

% Get base image
switch imageStr
    case 'Nucleus'
        imgname = [visualize.image_dir,filesep,eval(visualize.parameters.NucleusExpr)];
    case 'Cell'
        imgname = [visualize.image_dir,filesep,eval(visualize.parameters.CellExpr)];
    otherwise
        imgname = [visualize.image_dir,filesep,eval(visualize.parameters.(imageStr(1:end-3)).ImageExpr)];
end
image1 = checkread(imgname,visualize.parameters.BitDepth);
clim = prctile(image1(:),[1 99]);

% Get masks and overlay object (if it's in range).
load([visualize.label_dir,filesep,'xy',num2str(i),filesep,'NuclearLabels',filesep,'NuclearLabel-',numseq(frameNo,4),'.mat'])
image1(imdilate(NuclearLabel==cellno,ones(3))&~(NuclearLabel==cellno)) = min(image1(:));
cell_dir = [visualize.label_dir,filesep,'xy',num2str(i),filesep,'CellLabels'];

if exist(cell_dir,'dir')
    load([cell_dir,filesep,'CellLabel-',numseq(frameNo,4),'.mat'])
    image1(imdilate(CellLabel==cellno,ones(3))&~(CellLabel==cellno)) = max(image1(:));    
end
if size(image1,1)>size(image1,2)
   image1 = imrotate(image1,90);
end
imagesc(image1,'Parent',visualize.axes1)
axis image
colormap(gray)
set(visualize.axes1,'YTick',[],'XTick',[],'CLim',clim)
title(visualize.axes1,visualize.title)
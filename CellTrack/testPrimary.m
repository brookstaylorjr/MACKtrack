function [] = testPrimary(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% TESTPRIMARY Once parameters are loaded/set, segment specified images and display diagnostic figure
% This function is identical to TESTIMAGES, except specific to nuclear segmentation only.
%
% handles  master structure with parameters and naming data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% - - - - Setup: create structure for images, define options - - - -
p = handles.parameters;
locations = handles.locations;
i = min(p.XYRange);
j = min(p.TimeRange);

% Load all base images
imfo = imfinfo([locations.scope,p.ImagePath,eval(p.NucleusExpr)]);
bit_depth = imfo.BitDepth;
images = struct;
images.nucleus = checkread([locations.scope,p.ImagePath,eval(p.NucleusExpr)],bit_depth);

% Add measurement images
modules = p.ModuleNames;
for ind1 = 1:length(modules)
    % If the module is in use, and auxImageExpr is a file add it to the list of images
    if p.(modules{ind1}).Use 
        try
            if exist([locations.scope,p.ImagePath,filesep,eval(p.(modules{ind1}).ImageExpr)],'file')
                imgname = [locations.scope,p.ImagePath,filesep,eval(p.(modules{ind1}).ImageExpr)];
                images.(modules{ind1}) = checkread(imgname,bit_depth);
            end
        catch me
        end
    end
end
% Rotate images if necessary
imgnames = fieldnames(images);
for k = 1:length(imgnames)
    if size(images.(imgnames{k}),1)>size(images.(imgnames{k}),2)
        images.(imgnames{k}) = imrotate(images.(imgnames{k}),90);
    end
end

% - - - - 1st label: make label matrix from nuclear image - - - -
[data, diag_tmp] = primaryID(images.nucleus,p);
% Save information
handles.diagnostics.nuclei = diag_tmp;


assignin('base', 'diagnostics',handles.diagnostics)
assignin('base','images',images);
assignin('base','trackdata',data);


% - - - - Make all display images - - - - 
display_list = {};

% Make overlay list: nuclei, cells, and all primary measured images
overlayList = fieldnames(images);
% Create overlay list images
for k = 1:length(overlayList)
    % Get base image
    disp_img = modebalance(images.(overlayList{k}),1,bit_depth);
    x = [-2 5];
    disp_img(disp_img<x(1)) = x(1);
    disp_img(disp_img>x(2)) = x(2);
    disp_img = (disp_img - x(1))/(x(2)-x(1))*255;  
    % Overlay nuclear borders as orange
    border1 = (imdilate(data.nuclei,ones(3))-data.nuclei)>0;
    R = disp_img; R(border1) = R(border1)*0.25 + 0.75*248;
    G = disp_img; G(border1) = G(border1)*0.25 + 0.75*152;
    B = disp_img; B(border1) = B(border1)*0.25 + 0.75*29;
    disp_img = cat(3,R,G,B);
    disp_img = uint8(round(disp_img));
    handles.overlays.(overlayList{k}) = disp_img;
    display_list = cat(2,display_list,['overlay-',overlayList{k}]);
end

% - - - - Make full list of all diagnostic data - - - -
infoFields = fieldnames(handles.diagnostics);
for m = 1:length(infoFields)
   dataFields = fieldnames(handles.diagnostics.(infoFields{m}));
   for n = 1:length(dataFields)
       dataSize = size(handles.diagnostics.(infoFields{m}).(dataFields{n}));
       if (dataSize(1) == size(images.nucleus,1)) && (dataSize(2) == size(images.nucleus,2))
           display_list = cat(2,display_list,[infoFields{m},'-',dataFields{n}]);
       end
   end
end

handles.display_list = display_list;

% - - - - Initialize figure + popup - - - - 
% Initialize axes and slider and set properties
handles.diagnosticFig = figure;
handles.diagnosticAxes = axes('Parent',handles.diagnosticFig);
handles.diagnosticPopup  = uicontrol('Style', 'popupmenu','String',display_list,'Callback',{@popup_Callback,handles},'FontName','sans serif',...
    'BackgroundColor',[.99 .99 1]);
imagesc(handles.overlays.(overlayList{1}),'Parent',handles.diagnosticAxes), colormap gray
set(handles.diagnosticAxes,'YTick',[],'XTick',[])
set(handles.diagnosticFig,'ResizeFcn',{@fig_resize,handles},'Toolbar','figure');
assignin('base','p',handles.parameters)
guidata(handles.figure1,handles)


function popup_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Redraw graph if new variable is selected
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
newStr = handles.display_list{get(hObject,'Value')};
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
structType = newStr(1:strfind(newStr,'-')-1);
switch structType
    case 'overlay'
        imshow(handles.overlays.(newStr(strfind(newStr,'-')+1:end)),'Parent',handles.diagnosticAxes), colormap gray
    otherwise
        imagesc(handles.diagnostics.(structType).(newStr(strfind(newStr,'-')+1:end)),'Parent',handles.diagnosticAxes), colormap jet
end
axis image
set(gca,'xlim',xlim,'ylim',ylim)
set(handles.diagnosticAxes,'YTick',[],'XTick',[])

% ========================================================================================

function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize popup menu smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.diagnosticFig,'Position');
h = figPos(4);
w = figPos(3);
set(handles.diagnosticPopup,'Position',[floor(w/2)-200 5 400 24])
axis image
set(handles.diagnosticAxes,'OuterPosition',[0 30/h 1 (h-42)/h],'LooseInset',get(handles.diagnosticAxes,'TightInset')+[10/w 5/h 10/w 5/h])
% ========================================================================================
function [] = testImages(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = testImages(handles)
%
% TESTIMAGES Once parameters are loaded/set, segment specified images and display diagnostic figure
%
% handles  master structure with parameters and naming dat
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% - - - - Setup: create structure for images, define options - - - -
p = handles.parameters;
locations = handles.locations;
i = min(p.XYRange);
j = min(p.TimeRange);
p.i = i; p.j = j;

% Load all base images
imfo = imfinfo([locations.scope,p.ImagePath,eval(p.NucleusExpr)]);
bit_depth = imfo.BitDepth;

images = struct;
images.nucleus = checkread([locations.scope,p.ImagePath,eval(p.NucleusExpr)],bit_depth);
if ~strcmpi(p.ImageType,'none')
    images.cell = checkread([locations.scope,p.ImagePath,eval(p.CellExpr)],bit_depth);
else
    images.cell = images.nucleus;
end

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
% Rotate images (to make wider than tall)
imgnames = fieldnames(images);
rotate_flag = 0;
for k = 1:length(imgnames)
    if size(images.(imgnames{k}),1)>size(images.(imgnames{k}),2)
        images.(imgnames{k}) = imrotate(images.(imgnames{k}),90);
        rotate_flag = 1;
    end
end

% Convert any parameter flatfield images to functions
if isfield(p,'Flatfield')
    X = [];
    warning off MATLAB:nearlySingularMatrix
    for i = 1:length(p.Flatfield)
        if size(X,1) ~= numel(p.Flatfield{i})
            X = backgroundcalculate(size(p.Flatfield{i}));
        end        
        corr_img = p.Flatfield{i};
        pStar = (X'*X)\(X')*corr_img(:);
        % Apply correction
        corr_img = reshape(X*pStar,size(corr_img));
        if rotate_flag
            corr_img = imrotate(corr_img,90);
        end
        p.Flatfield{i} = corr_img-min(corr_img(:));
    end
end


switch lower(p.ImageType)
    case 'dic'
        % Set function names
        fnstem = 'dic';
        X = []; 
    case 'phase'
        fnstem = 'phase';
        X = backgroundcalculate(p.ImageSize);
    case 'none'
        fnstem = 'primary';
        X = [];
    case 'fluorescence'
        fnstem = 'fluorescence';
        X = images.nucleus; 
        
end

% Turn off combine structures warning
warning('off','MATLAB:combstrct')

% - - - - 1st mask: make mask from cell image - - - - 
tic
maskfn = str2func([fnstem,'ID']);
[data,diag_tmp] = maskfn(images.cell,p,X);
% Save information
handles.diagnostics.cell = diag_tmp;
tocs.CellMasking = toc;

% - - - - 1st label: make label matrix from nuclear image - - - -
tic
[data_tmp, diag_tmp] = nucleusID(images.nucleus,p,data);
% Save information
handles.diagnostics.nuclei = diag_tmp;
data = combinestructures(data_tmp,data);
tocs.NucMasking = toc;


% - - - - Check cells (misc check for binucleates, missed nuclei, etc.) - - - -
tic
[data_tmp,diag_tmp] = doubleCheck(data,images.cell,p);
% Save information
handles.diagnostics.check = diag_tmp;
data = combinestructures(data_tmp,data);
tocs.CheckCells = toc;

% - - - - 2nd label: make label matrix by segmenting cell image - - - - 
tic
segmentfn = str2func([fnstem,'Segment']);
[data_tmp, diag_tmp] = segmentfn(data,images.cell,p);
% Save information
handles.diagnostics.segmentation = diag_tmp;
data = combinestructures(data_tmp,data);
tocs.Segmentation = toc;

% Display times
n = fieldnames(tocs);
for k = 1:length(n)
    disp([n{k},'- ',num2str(tocs.(n{k})),' sec'])
end


assignin('base', 'diagnostics',handles.diagnostics)
assignin('base','images',images);
assignin('base','trackdata',data);


% - - - - Make all display images - - - - 
display_list = {};


setcolors;
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
    
    if ~strcmpi(p.ImageType,'none')
        %  Overlay cell borders as light blue
        border1 = (imdilate(data.cells,ones(3))-data.cells)>0;
        R(border1) = R(border1)*0.25 + 0.75*118;
        G(border1) = G(border1)*0.25 + 0.75*180;
        B(border1) = B(border1)*0.25 + 0.75*203;
    end
    
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
handles.imageSize = size(handles.overlays.(overlayList{1}));
handles.diagnosticFig = figure;
handles.diagnosticAxes = axes('Parent',handles.diagnosticFig);
handles.diagnosticPopup  = uicontrol('Style', 'popupmenu','String',display_list,'Callback',{@popup_Callback,handles},'FontName','sans serif',...
    'BackgroundColor',[.99 .99 1]);
% Make "Reset zoom" button
handles.reset  = uicontrol('Style', 'pushbutton','String','Reset zoom','BackgroundColor',[.99 .99 1]);
set(handles.reset,'Callback',{@reset_Callback,handles});
% Display image
imagesc(handles.overlays.(overlayList{1}),'Parent',handles.diagnosticAxes), colormap gray
set(handles.diagnosticAxes,'YTick',[],'XTick',[])
set(handles.diagnosticFig,'ResizeFcn',{@fig_resize,handles},'Toolbar','figure');
assignin('base','p',handles.parameters)
if isfield(handles,'figure1')
    guidata(handles.figure1,handles)
end
fig_resize([],[],handles)

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
set(handles.reset,'Position',[floor(w/2)+210,12, 70, 18]);
axis image
set(handles.diagnosticAxes,'OuterPosition',[0 30/h 1 (h-42)/h],'LooseInset',get(handles.diagnosticAxes,'TightInset')+[10/w 5/h 10/w 5/h])
% ========================================================================================


function reset_Callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Reset zoom state of figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

xlim = [1 handles.imageSize(2)];
ylim = [1 handles.imageSize(1)];

set(gca,'xlim',xlim,'ylim',ylim)

% ========================================================================================
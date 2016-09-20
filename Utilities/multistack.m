function [ ] = multistack(inputStack, bounds)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MULTISTACK Create a GUI that allows the user to see a 3-D stack of images using a slider at bottom of figure.
%
% inputStack     a 3-D matrix,  e.g. formed by cat(3,img1,img2...)
% bounds         (optional) max/min for image display
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin>1
    handles.flag1 = 1;
    handles.bounds = bounds;
else
    handles.flag1 = 0;
end


% Rotate image if it's taller than it is wide
if size(inputStack,1)>size(inputStack,2)
    inputStack = imrotate(inputStack,90);
end



clf(gcf,'reset')
zlength = size(inputStack,3);
handles.input = inputStack;

% Create the figure
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
% Display graph
set(handles.axes1,'YTick',[],'XTick',[])
set(handles.figure1,'ResizeFcn',{@fig_resize,handles},'Toolbar','figure');
slider_Callback(handles.slider1,[],handles)
fig_resize([],[],handles)


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
set(handles.text1,'String',[num2str(sliderVal),'/',num2str(size(handles.input,3))])
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
function [] = zeissbrowse(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = zeissbrowse(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ZEISSBROWSE creates GUI allowing user to browse a sequential image set inside a directory
% Image names MUST be formatted to specify channel, position, and timepoint as shown (e.g.):
% - fluorescence channel: c2 (must be single digit)
% - XY (stage) position: s05 (can be any number of digits)
% - timepoint: t001 (can be any number of digits)
%
% INPUT
% start_dir     image directory to start from (a dialog box will open if left unspecified)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<1
    if ispc; userdir= getenv('USERPROFILE'); 
    else userdir= getenv('HOME');
    end
    handles.dir = uigetdir(userdir,'Choose Image Folder');
else
    handles.dir = start_dir;
end

% Create the figure and axes
clf('reset')
handles.figure1 = gcf;
handles.axes1 = axes('Parent',handles.figure1,'YTick',[],'XTick',[]);
% Create status text
handles.text1 = uicontrol('Style','text','String','','BackgroundColor',[1 1 1],...
    'HorizontalAlignment','center');
handles.text2 = uicontrol('Style','text','String','','BackgroundColor',[1 1 1],...
    'HorizontalAlignment','left');

% Create slider, play button, and dropdown
handles.slider1  = uicontrol('Style','slider','BackgroundColor',[1 1 1]);
handles.dropdown1 = uicontrol('Style','popupmenu','String',{'Time (t)','Position (s)','Channel (c)'},...
    'BackgroundColor',[.99 .99 1]);
handles.histbutton = uicontrol('Style','pushbutton', 'BackgroundColor',[1 1 1],...
    'String', 'hist');
% Set resize and close functions
set(handles.figure1,'ResizeFcn',{@fig_resize,handles},'Toolbar','figure');
set(handles.figure1,'CloseRequestFcn',{@close_fig,handles},'Toolbar','figure');

% Get first image (jpg, tif) in the full list of files, and determine ranges of channel/position/time
ok_go = 0;
while ~ok_go
    try
        tic
        handles.dir_contents = cell(java.io.File(handles.dir).list());
        loadtoc = toc;
        ok_go = 1;
    catch ME
        handles.dir = uigetdir(start_dir,'Invalid folder - choose another one');
        if handles.dir==0
            error('No folder selected - exiting.')
        else
            tic
            handles.dir_contents = cell(java.io.File(handles.dir).list());
            loadtoc = toc;
            ok_go = 1;
        end
    end
end

disp(['Loaded directory contents in ',num2str(round(loadtoc*100)/100),' sec'])
ext_matches = {'.jpg', '.tif'};
draw = 0;
% Cycle through to get 1st image - parse its name, then break the loop
for i=1:length(handles.dir_contents)
    % First name found - parse name to find positions of time/xy/channel
    tmp_name = handles.dir_contents{i};
    ext_idx = find(strcmpi(tmp_name(max([1,end-3]):end), ext_matches),1);
    if ~isempty(ext_idx)
        handles.img = tmp_name;
        % Channel (c)
        [~,c1] = regexp(handles.img,'c[0-9]');
        c2 = regexp(handles.img(c1(end):end),'[^0-9]','once');
        handles.pos_c = c1(end):c1(end)+c2-2;
        handles.min_c = eval(handles.img(handles.pos_c));
        handles.max_c = handles.min_c;
        % Time (t)
        [~,t1] = regexp(handles.img,'t[0-9]');  
        t2 = regexp(handles.img(t1(end):end),'[^0-9]','once');
        handles.pos_t = t1(end):t1(end)+t2-2;
        handles.min_t = eval(handles.img(handles.pos_t));
        handles.max_t = handles.min_t;
        % Position (s)
        [~,s1] = regexp(handles.img,'s[0-9]');
        s2 = regexp(handles.img(s1(end):end),'[^0-9]','once');
        handles.pos_s = s1(end):s1(end)+s2-2;
        handles.min_s = eval(handles.img(handles.pos_s));
        handles.max_s = handles.min_s;       
        match_ind = 1:min([c1,t1,s1])-1;
        draw = 1;
        break
    end
end

if draw
    % Cycle the rest of the images and add to range (if necessary)
    for j=1:length(handles.dir_contents)           
        tmp_name = handles.dir_contents{j};
        match_ind(match_ind>length(tmp_name)) = [];
        if strcmp(tmp_name(match_ind), handles.img(match_ind))
            img = tmp_name;
            c_test = eval(img(handles.pos_c));
            if c_test > handles.max_c;
                handles.max_c = c_test;
            end
            if c_test < handles.min_c;
                handles.min_c = c_test;
            end
            % Time (t)
            t_test = eval(img(handles.pos_t));
            if t_test > handles.max_t;
                handles.max_t = t_test;
            end
            if t_test < handles.min_t;
                handles.min_t = t_test;
            end
            % Position (s)
            s_test = eval(img(handles.pos_s));
            if s_test > handles.max_s;
                handles.max_s = s_test;
            end
            if s_test < handles.min_s;
                handles.min_s = s_test;
            end                
        end
    end
    
    handles.img(handles.pos_s) = numseq(handles.min_s,length(handles.pos_s));
    handles.img(handles.pos_t) = numseq(handles.min_t,length(handles.pos_t));
    handles.img(handles.pos_c) = numseq(handles.min_c,length(handles.pos_c));

    
    % Set up callbacks
    set(handles.slider1,'Callback',{@gui_callback,handles});
    set(handles.dropdown1,'Callback',{@gui_callback,handles});
    set(handles.histbutton,'Callback',{@hist_callback,handles});

    % Set GUI elements to appropriate position
    set(handles.slider1,'Min',handles.min_t,'Max',handles.max_t, 'Value',eval(handles.img(handles.pos_t)), ...
        'SliderStep',[1/(handles.max_t-handles.min_t) 10/(handles.max_t-handles.min_t)]);
    set(handles.text1,'String', ['time:',handles.img(handles.pos_t),' | XY pos:', handles.img(handles.pos_s)...
        ' | channel:', handles.img(handles.pos_c)])
    set(handles.text2,'String',[numseq(eval(handles.img(handles.pos_t)),length(handles.pos_t)),'/',numseq(handles.max_t,length(handles.pos_t))])

    % Set starting image range (for display)
    imfo = imfinfo([handles.dir,filesep,handles.img]);
    handles.bit_depth = imfo.BitDepth;
    
    img = checkread([handles.dir,filesep,handles.img],handles.bit_depth,0,0);
    handles.CLim = double([min(img(:)), max(img(:))]);
    set(handles.axes1,'CLim',handles.CLim)
    if size(img,1) > size(img,2)
        img = imrotate(img,90);
    end
    image('CData',img,'Parent',handles.axes1,'CDataMapping','scaled','CData',img(end:-1:1,:)), colormap gray
    fig_resize([],[],handles)
    guidata(handles.figure1,handles)

else
    close(handles.figure1)
    error('Folder contains no valid images')
end

% ========================================================================================



function [handles_out] = drawimage(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Update image in the axes
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
img = checkread([handles.dir,filesep,handles.img],handles.bit_depth,0,0);
if size(img,1) > size(img,2)
    img = imrotate(img,90);
end
cla(handles.axes1)
image('Parent',handles.axes1,'CDataMapping','scaled','CData',img(end:-1:1,:))

if isfield(handles,'CLim')
    set(handles.axes1, 'CLim', handles.CLim)
else
    set(handles.axes1, 'CLim', [min(img(:)),max(img(:))])
    handles.CLim = double([min(img(:)),max(img(:))]);
end
colormap gray
handles_out = handles;
% ========================================================================================


function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize figure and  GUI elements smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.figure1,'Position');
h = figPos(4);
w = figPos(3);
set(handles.text1,'Position',[ceil(w/2)-200,h-24, 400,16]) % Centered at top
% Left-to-right bottom row elements
set(handles.histbutton,'Position',[4 4 50 22])
set(handles.slider1,'Position',[floor(w/2)-200 5 400 20])
set(handles.dropdown1,'Position',[floor(w/2)+210 9 100 18])
set(handles.text2,'Position', [floor(w/2)+320 6 120 16])
axis image
set(handles.axes1,'OuterPosition',[0 20/h 1 (h-40)/h],'LooseInset',get(handles.axes1,'TightInset')+[10/w 10/h 10/w 10/h])
set(handles.axes1,'YTick',[],'XTick',[])

% ========================================================================================

function gui_callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Callback for slider and dropdown menus - cycle position, channel, or timepoint
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Get dropdown selection: position, time, or channel
handles = guidata(handles.figure1);
list1 = get(handles.dropdown1,'String');
str1 = list1{get(handles.dropdown1,'Value')};
ltr = str1(end-1);
max1 = handles.(['max_',ltr']);
min1 = handles.(['min_',ltr']);
pos1 = handles.(['pos_',ltr]);
% Get/set slider value
if max1==get(handles.slider1,'Max')
    value1 = round(get(handles.slider1,'Value'));
else
    value1 = eval(handles.img(pos1));
end


% Set image name, slider and status texts
handles.img(pos1) = numseq(value1,length(pos1));
set(handles.slider1,'Max',max1,'Min',min1, 'Value',value1, 'SliderStep',[1/(max1-min1) 10/(max1-min1)])
set(handles.text1,'String', ['time:',handles.img(handles.pos_t),' | XY pos:', handles.img(handles.pos_s)...
    ' | channel:', handles.img(handles.pos_c)])
set(handles.text2,'String',[numseq(value1,length(pos1)),'/',numseq(max1,length(pos1))])

% If we're switching channel, redraw the histogram and reset CLim
if strcmp(ltr,'c')
    if isfield(handles,'CLim')
        handles = rmfield(handles,'CLim');
    end
    handles = drawhistogram(handles,1); % Redraw histogram    
end

% Load/draw image
handles = drawimage(handles);

guidata(handles.figure1, handles)
% ========================================================================================

function hist_callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Make ImageJ-styled figure to allow contrast adjustment
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles = guidata(handles.figure1);
if isempty(findobj('type','figure','name',['Adjust Contrast (', num2str(handles.figure1),')']))
    pos1 = get(handles.figure1,'Position');
    handles.hist = figure('Name', ['Adjust Contrast (', num2str(handles.figure1),')'],...
        'Position',[pos1(1)+pos1(3)+10 pos1(2)+pos1(4)-300 300 300], 'MenuBar','none', 'ToolBar','none');
end
handles.histaxes = axes('Parent',handles.hist);
cla(handles.histaxes,'reset')
% Get image bit depth and histogram
info1 = imfinfo([handles.dir,filesep,handles.img]);
handles.BitDepth = info1.BitDepth;
handles = drawhistogram(handles);
figPos = get(handles.hist,'Position');
h = figPos(4);
w = figPos(3);
set(handles.histaxes,'OuterPosition',[0 120/h 1 (h-120)/h],...
    'LooseInset',get(handles.histaxes,'TightInset')+[10/h 10/h 10/w 10/h])
% Create slider 1: image min
handles.histtext1 = uicontrol('Parent', handles.hist, 'Style','text','BackgroundColor',[1 1 1],'String','Min',...
    'Position',[10 100 w-10 16],'HorizontalAlignment','center');
handles.histslider1 = uicontrol('Parent', handles.hist, 'Style','slider','BackgroundColor',[1 1 1],...
    'Min',0,'Max',(2^handles.BitDepth)-1,'Value',handles.CLim(1),...
    'Position',[10 80 w-15 20],'Callback',{@histslider1_callback,handles});

% Create slider 2: image max
handles.histtext2 = uicontrol('Parent', handles.hist, 'Style','text','BackgroundColor',[1 1 1],'String','Max',...
    'Position',[10 60 w-10 16],'HorizontalAlignment','center');
handles.histslider2 = uicontrol('Parent', handles.hist, 'Style','slider','BackgroundColor',[1 1 1],...
    'Min',0,'Max',(2^handles.BitDepth)-1,'Value',handles.CLim(2),...
    'Position',[10 40 w-15 20],'Callback',{@histslider2_callback,handles});
guidata(handles.figure1, handles)
% ========================================================================================

function histslider1_callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Histogram slider 1: display minimum
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Get value/ update GUI elements
handles = guidata(handles.figure1);
newVal = round(get(hObject,'Value'));
if newVal>=handles.CLim(2)
    newVal = handles.CLim(2)-1;
end
handles.CLim(1) = newVal;
drawnow
set(handles.axes1,'CLim',handles.CLim)
set(hObject,'Value',newVal)
% Redo graph annotation
hold(handles.histaxes,'on')
if isfield(handles,'histline1')
    delete(handles.histline1)
end
handles.histline1 = plot(handles.histaxes,ones(1,2)*handles.CLim(1),get(handles.histaxes,'YLim'),...
    'LineWidth', 3, 'Color',[.8 .5 .5]);
hold(handles.histaxes,'off')
guidata(handles.figure1,handles)
% ========================================================================================

function histslider2_callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Histogram slider 2: display maximum
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Get value/ update GUI elements
handles = guidata(handles.figure1);
newVal = round(get(hObject,'Value'));
if newVal<=handles.CLim(1)
    newVal = handles.CLim(1)+1;
end
handles.CLim(2) = newVal;
drawnow
set(handles.axes1,'CLim',handles.CLim)
set(hObject,'Value',newVal)
% Redo graph annotation
hold(handles.histaxes,'on')
if isfield(handles,'histline1')
    delete(handles.histline2)
end
handles.histline2 = plot(handles.histaxes,ones(1,2)*handles.CLim(2),get(handles.histaxes,'YLim'),...
    'LineWidth', 3, 'Color',[.8 .5 .5]);
hold(handles.histaxes,'off')
guidata(handles.figure1,handles)
% ========================================================================================


function handles_out = drawhistogram(handles,flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Draw the actual relative frequency image histogram
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ~isempty(findobj('type','figure','name',['Adjust Contrast (', num2str(handles.figure1),')']))
    cla(handles.histaxes,'reset')
    if nargin<2
        flag = 0;
    end
    img = double(imread([handles.dir,filesep,handles.img],'TIFF'));
    x = 0:2^(handles.BitDepth-8):(2^handles.BitDepth)-1;
    n = hist(img(:), x);
    n = n/sum(n);
    bar(handles.histaxes,x,n, 1, 'LineStyle','none')
    axis(handles.histaxes,[min(x) max(x) 0 max(n)*1.1])
    if flag
        handles.CLim = [min(img(:)), max(img(:))];
    end
    
    % Show lines showing with existing thresholds
    hold(handles.histaxes,'on')
    handles.histline1 = plot(handles.histaxes,ones(1,2)*handles.CLim(1),get(handles.histaxes,'YLim'),...
        'LineWidth', 3, 'Color',[.8 .5 .5]);
    handles.histline2 = plot(handles.histaxes,ones(1,2)*handles.CLim(2),get(handles.histaxes,'YLim'),...
        'LineWidth', 3, 'Color',[.8 .5 .5]);
    hold(handles.histaxes,'off')

    % Set hist sliders
    if isfield(handles,'histslider1')
        set(handles.histslider1,'Value',handles.CLim(1))
        set(handles.histslider2,'Value',handles.CLim(2))
    end
end

handles_out = handles;
% ========================================================================================


function [] = close_fig(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Close histogram with figure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
delete(findobj('type','figure','name',['Adjust Contrast (', num2str(handles.figure1),')']))


if isempty(gcbf)
   if length(dbstack) == 1
      warning('MATLAB:closereq',...
      'Calling closereq from the command line is now obsolete, use close instead');
   end
   close force
else
   delete(gcbf);
end
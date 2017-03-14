function [] = browse3i(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = browse3i(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BROWSE3I creates GUI allowing user to browse a sequential image set inside a directory
% Image names MUST be formatted to specify channel, position, and timepoint as shown (e.g.):
% - fluorescence channel: "_C2" (must be single digit)
% - XY (stage) position: "Position 05" (can be any number of digits). NOT ZERO PADDED.
% - timepoint: "_T001" (can be any number of digits)
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
handles.dropdown1 = uicontrol('Style','popupmenu','String',{'Position (s)','Time (t)','Channel (c)'},...
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
        handles.dir_contents = quickdir(handles.dir);
        loadtoc = toc;
        ok_go = 1;
    catch ME
        handles.dir = uigetdir(start_dir,'Invalid folder - choose another one');
        if handles.dir==0
            error('No folder selected - exiting.')
        else
            tic
            handles.dir_contents = quickdir(handles.dir);
            loadtoc = toc;
            ok_go = 1;
        end
    end
end
disp(['Loaded directory contents in ',num2str(round(loadtoc*100)/100),' sec'])


% Need some kind of way to gracefully handle "extra" data points -> for now, filter them out.
filter_str = @(str) (length(str)>10) && strcmp(str(end-10),'.');
drop_names  = cellfun(filter_str,handles.dir_contents);
if sum(drop_names)>0
    disp('Warning: found "intermediate" timepoint images in this folder (e.g. ''t001.02_C0.tiff'') - dropping these.')
    handles.dir_contents(drop_names) = [];
end


ext_matches = {'.jpg', '.tif', 'tiff', '.png'};
draw = 0;
% Cycle through to get 1st image - parse its name, then break the loop
for i=1:length(handles.dir_contents)
    % First name found - parse name to find positions of time/xy/channel
    tmp_name = handles.dir_contents{i};
    ext_idx = find(strcmpi(tmp_name(max([1,end-3]):end), ext_matches),1);
    disp(['''',ext_matches{ext_idx}, ''' image found - scanning for more of these files...'])
    if ~isempty(ext_idx)
        handles.img = tmp_name;
        
        % String part 1: before XY position (s)
        [~,s1] = regexp(handles.img,'Position [0-9]');
        s2 = regexp(handles.img(s1(end):end),'[^0-9]','once');
        pos_s = s1(end):s1(end)+s2-2;
        end1 = min(pos_s)-1; 
        start2 = max(pos_s)+1;
        handles.string{1} = handles.img(1:end1);
        handles.min_s = eval(handles.img(pos_s));
        handles.max_s = handles.min_s;       
        
        % Time (t)
        [~,t1] = regexp(handles.img,'_T[0-9]');  
        t2 = regexp(handles.img(t1(end):end),'[^0-9]','once');
        pos_t = t1(end):t1(end)+t2-2;
        end2  = min(pos_t)-1; 
        start3 = max(pos_t)+1;
        handles.string{2} = handles.img(start2:end2);
        handles.min_t = eval(handles.img(pos_t));
        handles.max_t = handles.min_t;
        
        
        % Channel (c)
        [~,c1] = regexp(handles.img,'_C[0-9]');
        c2 = regexp(handles.img(c1(end):end),'[^0-9]','once');
        pos_c = c1(end):c1(end)+c2-2;
        end3 = min(pos_c)-1;
        start4 = max(pos_c)+1;
        handles.string{3} = handles.img(start3:end3);
        handles.string{4} = handles.img(start4:end);
        handles.min_c = eval(handles.img(pos_c));
        handles.max_c = handles.min_c;        
        
        handles.lengths = [length(pos_s), length(pos_t), length(pos_c)];
        draw = 1;
        break
    end
end


if draw
    % Cycle the rest of the images and add to range (if necessary)    
    for j=1:length(handles.dir_contents)           
        tmp_name = handles.dir_contents{j};
        ext_idx = find(strcmpi(tmp_name(max([1,end-3]):end), ext_matches),1);
        if ~isempty(ext_idx)
            img = tmp_name;
            [test_vals] = extract_vals(img(length(handles.string{1}):end)); 
            if test_vals(2) > handles.max_t;
                handles.max_t = test_vals(2);
            end
            if test_vals(2) < handles.min_t;
                handles.min_t = test_vals(2);
            end
            if test_vals(1) > handles.max_s;
                handles.max_s = test_vals(1);
            end
            if test_vals(1) < handles.min_s;
                handles.min_s = test_vals(1);
            end
            if test_vals(3) > handles.max_c;
                handles.max_c = test_vals(3);
            end
            if test_vals(3) < handles.min_c;
                handles.min_c = test_vals(3);
            end
        end
    end
    
    handles.vals = [handles.min_s, handles.min_t, handles.min_c];
    handles.img = form_img(handles.string, handles.vals, handles.lengths);
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
    
    % Set up callbacks
    set(handles.slider1,'Callback',{@gui_callback,handles});
    set(handles.dropdown1,'Callback',{@gui_callback2,handles});
    set(handles.histbutton,'Callback',{@hist_callback,handles});

    % Set GUI elements to appropriate position
    set(handles.slider1,'Min',handles.min_s,'Max',handles.max_s, 'Value',handles.vals(1), ...
        'SliderStep',[1/(handles.max_s-handles.min_s) 10/(handles.max_s-handles.min_s)]);
    set(handles.text1,'String', ['time:',num2str(handles.vals(2)),' | XY pos:', num2str(handles.vals(1))...
        ' | channel:', num2str(handles.vals(3))])
    set(handles.text2,'String',[numseq(handles.vals(1),handles.lengths(1)),'/',numseq(handles.max_s,handles.lengths(1))])
    
    
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
    set(handles.axes1, 'CLim', [min(img(:)),max(img(:))+1])
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
% Callback for slider - cycle position, channel, or timepoint
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Get dropdown selection: position, time, or channel
handles = guidata(handles.figure1);
list1 = get(handles.dropdown1,'String');
idx1 = get(handles.dropdown1,'Value');
str1 = list1{idx1};
ltr = str1(end-1);
max1 = handles.(['max_',ltr']);
min1 = handles.(['min_',ltr']);
% Get/set slider value

handles.vals(idx1) = round(get(handles.slider1,'Value'));

% Set image name, slider and status texts
handles.img = form_img(handles.string, handles.vals, handles.lengths);
set(handles.slider1,'Max',max1,'Min',min1, 'Value',handles.vals(idx1), 'SliderStep',[1/(max1-min1) 10/(max1-min1)])
set(handles.text1,'String', ['time:',num2str(handles.vals(2)),' | XY pos:', num2str(handles.vals(1))...
    ' | channel:', num2str(handles.vals(3))])
set(handles.text2,'String',[numseq(handles.vals(idx1),handles.lengths(idx1)),'/',numseq(max1,handles.lengths(idx1))])

% Load/draw image
handles = drawimage(handles);

guidata(handles.figure1, handles)

% If we're switching channel, redraw the histogram and reset CLim
if strcmp(ltr,'c')
    if isfield(handles,'CLim')
        handles = rmfield(handles,'CLim');
    end
    handles = drawhistogram(handles,1); % Redraw histogram    
end

guidata(handles.figure1, handles)
% ========================================================================================

function gui_callback2(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Callback for dropdown menu - cycle position, channel, or timepoint
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Get dropdown selection: position, time, or channel
handles = guidata(handles.figure1);
list1 = get(handles.dropdown1,'String');
idx1 = get(handles.dropdown1,'Value');
str1 = list1{idx1};
ltr = str1(end-1);
max1 = handles.(['max_',ltr']);
min1 = handles.(['min_',ltr']);

% Set image name, slider and status texts
handles.img = form_img(handles.string, handles.vals, handles.lengths);
set(handles.slider1,'Max',max1,'Min',min1, 'Value',handles.vals(idx1), 'SliderStep',[1/(max1-min1) 10/(max1-min1)])
set(handles.text1,'String', ['time:',num2str(handles.vals(2)),' | XY pos:', num2str(handles.vals(1))...
    ' | channel:', handles.vals(3)])
set(handles.text2,'String',[numseq(handles.vals(idx1),handles.lengths(idx1)),'/',numseq(max1,handles.lengths(idx1))])

guidata(handles.figure1, handles)
% ========================================================================================


function hist_callback(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Make ImageJ-styled figure to allow contrast adjustment
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
handles = guidata(handles.figure1);
if isfield(handles,'histaxes')
    handles = rmfield(handles,'hist');
    handles = rmfield(handles,'histaxes');
    handles = rmfield(handles,'histline1');
    handles = rmfield(handles,'histline2');
    handles = rmfield(handles,'histtext1');
    handles = rmfield(handles,'histtext2');
    handles = rmfield(handles,'histslider1');
    handles = rmfield(handles,'histslider2');
end
if isempty(findobj('type','figure','name',['Adjust Contrast (', num2str(double(handles.figure1)),')']))
    pos1 = get(handles.figure1,'Position');
    handles.hist = figure('Name', ['Adjust Contrast (', num2str(double(handles.figure1)),')'],...
        'Position',[pos1(1)+pos1(3)+10 pos1(2)+pos1(4)-300 300 300], 'MenuBar','none', 'ToolBar','none');
end
handles.histaxes = axes('Parent',handles.hist);
cla(handles.histaxes,'reset')
% Get image bit depth and histogram
info1 = imfinfo([handles.dir,filesep,handles.img]);
handles.BitDepth = info1.BitDepth;
guidata(handles.figure1, handles)
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
handles = guidata(handles.figure1);

if ~isempty(findobj('type','figure','name',['Adjust Contrast (', num2str(double(handles.figure1)),')']))
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
delete(findobj('type','figure','name',['Adjust Contrast (', num2str(double(handles.figure1)),')']))


if isempty(gcbf)
   if length(dbstack) == 1
      warning('MATLAB:closereq',...
      'Calling closereq from the command line is now obsolete, use close instead');
   end
   close force
else
   delete(gcbf);
end


function string_out = form_img(string_group, vals, lengths)
string_out = [string_group{1}, num2str(vals(1)), string_group{2}, numseq(vals(2),lengths(2))....
string_group{3}, numseq(vals(3),lengths(3)), string_group{4}];


function [vals] = extract_vals(string_in)
all_numidx = regexp(string_in,'[0-9]');
a = find(diff([2, all_numidx])~=1);
b = [a(2:end)-1 length(all_numidx)];
vals = [eval(string_in(all_numidx(a(1):b(1)))),...
    eval(string_in(all_numidx(a(2):b(2)))),...
    eval(string_in(all_numidx(a(3):b(3))))];
function handlesOut = setOptions(handles)

% Check for existence of preexisting options- set otherwise
if ~isfield(handles,'Options')
    handles.Options.DataDirectory = handles.locations.data;
    handles.Options.MinLifetime = 0;
    handles.Options.FilterFlags = [0 0 0 0];
    handles.Options.Smoothing = ' None';
    handles.Options.SmoothingWindow = 3;
    handles.Options.PixelConversion = 4.03; % pixels per micron, was 3.3 on 20x Nikon Te/Ti
    handles.Options.FramesPerHour = 12;
    handles.Options.GraphMin1 = nan;
    handles.Options.GraphMin2 = nan;
    handles.Options.GraphMax1 = nan;
    handles.Options.GraphMax2 = nan;
    handles.Options.LogCompress1 = 0;
    handles.Options.LogCompress2 = 0;
    handles.Options.Grouping = ' None';
    handles.Options.GroupSubpopulations = 0;
    handles.Options.GraphType = ' Line plot';
    handles.Options.Subpopulations = 2;
    handles.Options.Measure2Grouping = 0;
    handles.Options.HistogramBins1 = 30;
    handles.Options.HistogramBins2 = 30;
    handles.Options.HistogramMax1 = 0.5;
    handles.Options.HistogramMax2 = 0.5;
    handles.Options.StartTime = 0.0;
    handles.Options.VisualizeXY = 1;
    handles.Options.VisualizeCell  = 1;
end


checkbox_names = 'ABCDEFG';

% Set all GUI elements to reflect option set.

% uipanel1: Parent Directory
set(handles.edit1A,'String',handles.Options.DataDirectory)

% uipanel3: Filtering Criteria
for i = 1:4
        set(handles.(['checkbox3',checkbox_names(i)]),'Value',handles.Options.FilterFlags(i))
end
set(handles.edit3A,'String',num2str(handles.Options.MinLifetime))

% uipanel6: Data Smoothing
i = find(strcmp(handles.Options.Smoothing,{' None', ' Median',' Mean', ' Loess'}));
set(handles.uipanel6,'SelectedObject',handles.(['radiobutton6',checkbox_names(i)]))
set(handles.edit6A,'String',num2str(handles.Options.SmoothingWindow))

% uipanel7: Graph Options
set(handles.edit7A,'String',num2str(handles.Options.PixelConversion))
set(handles.edit7B,'String',num2str(handles.Options.FramesPerHour))
set(handles.edit7C,'String',num2str(handles.Options.GraphMin1))
set(handles.edit7D,'String',num2str(handles.Options.GraphMin2))
set(handles.edit7E,'String',num2str(handles.Options.GraphMax1))
set(handles.edit7F,'String',num2str(handles.Options.GraphMax2))
set(handles.checkbox7A,'Value',handles.Options.LogCompress1)
set(handles.checkbox7A,'Value',handles.Options.LogCompress2)
set(handles.edit7G,'String',num2str(handles.Options.HistogramBins1))
set(handles.edit7H,'String',num2str(handles.Options.HistogramBins2))
set(handles.edit7I,'String',num2str(handles.Options.HistogramMax2))
set(handles.edit7J,'String',num2str(handles.Options.HistogramMax2))
set(handles.edit7K,'String',numseq(fix(handles.Options.StartTime),2))
set(handles.edit7K_2,'String',numseq(round((handles.Options.StartTime-fix(handles.Options.StartTime))*60),2))


% uipanel8: Grouping type
i = strcmp(handles.Options.Grouping,{' None', ' By lineage',' Gaussian mix (1-D)'});
set(handles.uipanel8,'SelectedObject',handles.(['radiobutton8',checkbox_names(i)]))
set(handles.edit8A,'String', num2str(handles.Options.Subpopulations))
set(handles.checkbox8A,'Value',handles.Options.Measure2Grouping)


% uipanel9: Graph type
i = strcmp(handles.Options.GraphType,{' Line plot', ' Stack colormap',' Histogram series',' Scatter plot'});
set(handles.uipanel9,'SelectedObject',handles.(['radiobutton9',checkbox_names(i)]))


handlesOut = handles;
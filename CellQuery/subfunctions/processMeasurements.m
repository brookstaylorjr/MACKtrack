function handlesOut  = processMeasurements(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PROCESSCELLS  Pulls selected measurements, drops cells, smooths data, models (if selected),
% and sets graph options. ALL pulled data is saved under handles.Export
%
% handles      master structure with measurement data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 11/15/2011: initial commit
% 02/04/2012: separated from plotting command (so it can be folded into grouping)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Columns in CellData: 
% 1) xyPosition
% 2) index in xy pos
% 3) FrameIn 
% 4) FrameOut 
% 5) Parent 
% 6) EdgeCell
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Pull existing export data, add CellData and filter dropped cells
export = handles.Export;
export.CellData = handles.Measurements.CellData;
export.CellData(handles.DropCells,:) = [];
export.Measurements = struct;

% Process all fields (filter/smooth/convert to microns)
names = fieldnames(handles.Measurements);
disp(['Smoothing:',handles.Options.Smoothing,'...'])
for n = 1:length(names)
    if ~strcmp(names{n},'CellData') 
        measurement = handles.Measurements.(names{n});
        measurement(handles.DropCells,:) = [];
        if ~iscell(handles.Measurements.(names{n}))
            measurement = smoothMeasurement(measurement,handles.Options,export.CellData, names{n});
        end
        export.Measurements.(names{n}) = measurement;
    end
end


% Get graphed measurements 
for m = 1:2 
    % Create Measurement and CellData fields/names 
    field = ['MeasureField',num2str(m)];
    measurement = ['Measurement',num2str(m)];
    letters = {'A','B'};
    export.(measurement) =  export.Measurements.(handles.Export.(field));
    
    % Log compress data, if selected, and set graph limits  
    % Check data for log compression
    if ((handles.Options.(['LogCompress',num2str(m)])) && (min(export.(measurement)(:))<=0))
        disp(['Data out of range for measurement ', num2str(m),'; ignoring log compression'])
        handles.Options.(['LogCompress',num2str(m)]) = 0;
        set(handles.(['checkbox7',letters{m}]),'Value',0)
    end

    % Log-compressed case: use traditional logarithmic ticks
    if handles.Options.(['LogCompress',num2str(m)]) 
        export.(measurement) = log10(export.(measurement));
        export.([measurement,'Info']).LogCompress = 1;
    else
        export.([measurement,'Info']).LogCompress = 0;
    end
    
end

% Update Export field and perform grouping/statistical modeling on data

handles.Export = export;
if ~handles.Flags.Modeled
    handles = modelTrajectories(handles);
end

% Rearrange all measurements/celldata
handles.Export.CellData = handles.Export.CellData(handles.Export.Order,:);
handles.Export.Measurement1 = handles.Export.Measurement1(handles.Export.Order,:);
handles.Export.Measurement2 = handles.Export.Measurement2(handles.Export.Order,:);

% Filter, convert, and rearrange all measurements
for n = 1:length(names)
    if ~strcmp(names{n},'CellData')
        handles.Export.Measurements.(names{n}) = handles.Export.Measurements.(names{n})(handles.Export.Order,:);
    end
end

handlesOut = handles;
% ========================================================================================



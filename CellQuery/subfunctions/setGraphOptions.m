function handlesOut  = setGraphOptions(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SETGRAPHOPTIONS  Sets graph options. ALL  data is saved under handles.Export.Measurement(n)Info
%
% handles      master structure with measurement data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 11/15/2011: initial writing
% 02/04/2012: separated from plotting command (so it can be folded into grouping)

export = handles.Export;
for m = 1:2
    % Create Measurement and MeasurementInfo fields/names
    field = ['MeasureField',num2str(m)];
    measurement = ['Measurement',num2str(m)];
    info = ['Measurement',num2str(m),'Info'];
    export.(info) = struct;
    
    % Get time vector
    times = (0:size(export.(measurement),2)-1)/handles.Options.FramesPerHour - handles.Options.StartTime;

    % Get measurement bounds - use user-set value, or auto-set to data max/min
    log_comp =  handles.Options.(['LogCompress',num2str(m)]);
    
    if isnan(handles.Options.(['GraphMax',num2str(m)])) % No maximum set, get from data
        m_tmp = export.(measurement)(:);
        m_tmp(isinf(m_tmp)|isnan(m_tmp)) = [];
        max_val = 1.05*max(m_tmp);
    else
        if log_comp
            max_val = log10(handles.Options.(['GraphMax',num2str(m)]));
        else
            max_val = handles.Options.(['GraphMax',num2str(m)]);
        end    
    end
    
    if isnan(handles.Options.(['GraphMin',num2str(m)])) % No minimum set, get from data
        m_tmp = export.(measurement)(:);
        m_tmp(isinf(m_tmp)|isnan(m_tmp)) = [];
        min_val = 0.95*min(m_tmp);
    else
        if log_comp
            min_val = log10(handles.Options.(['GraphMin',num2str(m)]));
            if min_val<=0
                m_tmp = export.(measurement)(:);
                m_tmp(isinf(m_tmp)|isnan(m_tmp)) = [];
                min_val = 0.95*min(m_tmp);
            end
        else
            min_val = handles.Options.(['GraphMin',num2str(m)]);
        end
    end
    bounds = [min(min_val,max_val),max(min_val,max_val)];
    
    % Set up graph options - ticks, log compression, etc.
    export.(info) = maketicks(times,bounds,log_comp);
    
    % Add measurement name, unit, and log compression
    switch export.(field)
        case 'Area'
            units = ' (um sq.)';
        case {'Movement_Cell', 'Movement_Nucleus', 'Perimeter'} 
            units = ' (um)';
        otherwise
            units = '';
    end
    export.(info).Name = [export.(field), units];
    
end
handles.Export = export;
handlesOut = handles;
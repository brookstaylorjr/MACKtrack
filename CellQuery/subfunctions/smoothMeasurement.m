function measure_out = smoothMeasurement(measure_in, options, CellData, measure_name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SMOOTHMEASUREMENT filters, smooths, and converts pixel-based measurements to microns
%
% measure_in     input measurement array
% options        specified by CellQuery GUI - need Smoothing, SmoothingWindow, PixelConversion, FramesPerHour
% CellData       information on cell trajectories generated during tracking
% measure_name   measurement's name
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Loop rows: assign NaN before frame in/after frame out, smooth
for i = 1:size(measure_in,1) 
    if CellData(i,4)<size(measure_in,2)
        inFrames = CellData(i,3):(CellData(i,4)-round(options.FramesPerHour*0.75)); % Omit last 45min of frames
    else
        inFrames = CellData(i,3):CellData(i,4);
    end
    % Convert out-of-range frames to NaN
    nanFrames = 1:size(measure_in,2);
    nanFrames(inFrames)=[];
    measure_in(i,nanFrames) = nan;    
    % Smoothing
    switch options.Smoothing      
    case ' Median'
        measure_in(i,inFrames) = medfilt1(measure_in(i,inFrames),options.SmoothingWindow);
    case ' Mean'
        measure_in(i,inFrames) = smooth(measure_in(i,inFrames),options.SmoothingWindow,'moving');
    case ' Loess'
        measure_in(i,inFrames) = smooth(measure_in(i,inFrames),options.SmoothingWindow/length(inFrames),'loess');
    end 
end

% Divide measurements by calibration
switch measure_name
    case 'Area'
        measure_in = measure_in/(options.PixelConversion.^2);
    case {'Movement_Cell', 'Movement_Nucleus', 'Perimeter'} 
        measure_in = measure_in/(options.PixelConversion);
end
measure_out = measure_in;
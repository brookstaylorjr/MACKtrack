function [info] = maketicks(times, bounds, log_comp)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MAKETICKS(times, bounds, log_comp) creates tick locations/labels from a time/measurement
%
% times          vector of time points
% bounds         user-set graph max and minimum
% log_comp       assume log compression for measurement (1 or 0)
%
% info           output structure with all ticks and ticklabels
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Default is non-log compression
if nargin < 3
    log_comp = 0;
end

% Time ticks: try largest power of ten, or use integral step that leads to best fit with 4-7 ticks total

test_val = max(floor(times)) / 10^(floor(log10(max(floor(times)))));

if (test_val>3) && (test_val<7)
    info.TimeTicks = 0:10^(floor(log10(max(floor(times))))):max(floor(times));
else
    tick_nos = 3:6;
    a = max(floor(times))./tick_nos;
    b = find((a-floor(a))==min(a-floor(a)),1,'last');
    info.TimeTicks = 0:floor(a(b)):floor(a(b))*tick_nos(b);
end

info.TimeBounds = [min(times), max(times)];
info.Times = times;
info.LogCompress = log_comp;

% Measurement bounds
info.MeasurementBounds = bounds;

if log_comp % Log-compressed case: use traditional logarithmic ticks  
    log_min = min(bounds);
    log_max = max(bounds);
    % Find nearest power of ten above and below data to mark lower limit, mark other powers of ten
    % Make unlabeled tick marks for e.g. 20, 30... and 200, 300... etc
    i = 1;
    j = 0;
    ticks = floor(log_min);
    tick_names = {num2str(10^ticks(end))};

    while ticks(end)<log_max
        i = mod(i+1,10);
        if i==0
            j = j+1;         
        else
            ticks = [ticks ticks(1)+log10(i)+j];
            if i == 1
                tick_names = [tick_names;num2str(10^ticks(end))];
            else
                tick_names = [tick_names; ' '];
            end
        end
    end    
    info.MeasurementTicks = ticks;
    info.MeasurementTickLabels = tick_names;
   
% Non-log case: choose integral ticks, after normalizing range to a value from 5-50.
else 
    % Convert measument range to value from 5-50
    range_orig = max(bounds) - min(bounds);
    range_conv = range_orig/10^(floor(log10(range_orig)));
    if (range_conv<5) 
        range_conv = range_conv*10; 
    end
    range_factor = range_orig/range_conv;
    % Move up bottom tick, recalculate range
    new_min = ceil(min(bounds/range_factor))*range_factor;
    new_range = max(bounds) - new_min;
    % Optimize tick placement
    tick_nos = 4:7;
    a = (new_range/range_factor)./tick_nos;
    b = find((a-floor(a))==min(a-floor(a)),1,'last');
    info.MeasurementTicks = new_min + (0:range_factor*floor(a(b)):range_factor*floor(a(b))*tick_nos(b));
    info.MeasurementTickLabels = cellstr(num2str(info.MeasurementTicks'));
end

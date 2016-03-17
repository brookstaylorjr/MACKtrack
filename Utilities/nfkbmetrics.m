function [metrics,aux, graph, info, measure] = nfkbmetrics(id,varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% metrics = nfkbmetrics(id)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMETRICS uses the see_nfkb_native function to filter and preprocess NFkB trajectories,
% then calculates related metrics regarding activation. Metric types include:
% 
% 1) time series (base NFkB dynamics, resampled to 12 frames/hr
% 2) integrated activity
% 3) differentiated activity
% 4) calculated metrics: measuring aspects of oscillation, duration, timing ,and amplitude
%
% INPUTS (required):
% id             filename or experiment ID (from Google Spreadsheet specified in "locations.mat")
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Display'      'on' or 'off' - show graphs (default: process data only; no graphs)
% 'Verbose'      'on' or 'off' - show verbose output
% 'Endframe'     final frame used to filter for long-lived cells (default = 100)
%
% OUTPUT: 
% metrics   structure with output fields
% aux       Extra data (e.g. fourier information (FFT, power, frequencies), thresholds used in envelope/duration)
% graph     main structure output from see_nfkb_native
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
% Required: ID input
valid_id = @(x) assert((isnumeric(x)&&length(x)==1)||exist(x,'file'),...
    'ID input must be spreadsheet ID or full file path');
addRequired(p,'id',valid_id);
% Optional parameters
addParameter(p,'Baseline', 1.9,@isnumeric);
addParameter(p,'Endframe',100, @isnumeric);
parse(p,id, varargin{:})

%% PARAMETETERS for finding off times - chosen using 'scan_off_params.m'
baseline = p.Results.Baseline; % Minimum activity required for cell to register as 'on'
window_sz = 14; % 1+ hr windows (on either side of a given timepoint)
thresh = 0.9; % Pct of inactivity allowed in a given window
cutoff_time = 4; % time to look for cell activity before declaring it "off" (hrs)
off_pad = 6; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)

%% INITIALIZATION. Load and process data. Interpolate time series, calculate deriv/integral approximations
if ismember('Endframe',p.UsingDefaults)
    [graph, info, measure] = see_nfkb_native(id);
else
   [graph, info, measure] = see_nfkb_native(id,'Endframe',p.Results.Endframe);
   graph.var = graph.var(:,1:p.Results.Endframe);
   graph.t = graph.t(1:size(graph.var,2));
end

% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
t = min(graph.t):1/12:max(graph.t);
if length(t)~=length(graph.t)
    metrics.time_series = nan(size(graph.var,1),length(t));
    for i = 1:size(graph.var,1)
        metrics.time_series(i,:) = interp1(graph.t,graph.var(i,:),t);
    end
else
    metrics.time_series = graph.var;
end

% 2) integrated activity
metrics.integrals = nan(size(metrics.time_series));
for i = 1:size(metrics.integrals,1)
    metrics.integrals(i,:) = cumtrapz(t,metrics.time_series(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = medfilt1(metrics.time_series,3,[],1);
metrics.derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);

% 4) calculated metrics
% MAX/MIN metrics
metrics.max_amplitude = nanmax(metrics.time_series,[],2);
metrics.max_integral = nanmax(metrics.integrals,[],2);
metrics.max_derivative = nanmax(metrics.derivatives,[],2);
metrics.min_derivative = nanmin(metrics.derivatives,[],2);

%% DURATION(1) Compute an off-time for all cells
metrics.off_times = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];
inactive = smoothrows(inactive<(baseline),(window_sz*2));
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1;
inactive = inactive(:,frontcrop:end);
inactive = inactive(:,1:size(smoothed,2));
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*12) % ignore cells who only turned on after 6+ hrs.
            metrics.off_times(i) = active_times(end);
        end
    end    
end
metrics.off_times = (metrics.off_times-1)/12;
metrics.off_times(metrics.off_times<0) = 0;

%% METRICS OF OSCILLATION
% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times)*12;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series,1)
    if(metrics.off_times(i)>0)
        y = metrics.time_series(i,1:depth);
        off_frame = min([length(y), metrics.off_times(i)*12+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
%     % (Code to check this "second-harmonic" thing)
%     if i<49
%     figure('Position',positionfig(220,100,[6,3])),
%     ha = tight_subplot(1,2);
%     plot(ha(1),1:100,metrics.time_series(i,1:100))
%     set(ha(1),'Ylim',[0 9],'XLim',[0 100],'Box','on')
%     hold(ha(2),'on')
%     plot(ha(2),freq, aux.power(i,:))
%     plot(ha(2), freq(locs),pks,'o')
%     hold(ha(2),'off')
%     set(ha(2),'XLim',[0 2],'Box','on')
%     end
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq(i) = 3600*aux.freq(1);
    end
end
%%
% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series,1)
        metrics.oscfrac(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac(i,j))
            metrics.oscfrac(i,j) = 0;
        end
    end
end



%% METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude
metrics.pk1_time = nan(size(metrics.time_series,1),1);
metrics.pk1_amp =  nan(size(metrics.time_series,1),1);
metrics.pk2_time = nan(size(metrics.time_series,1),1);
metrics.pk2_amp =  nan(size(metrics.time_series,1),1);
for i = 1:size(metrics.pk1_time,1)    
    [pks, locs] = globalpeaks(metrics.time_series(i,1:min([90,p.Results.Endframe])),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];  
    end
    pks(locs<4) = [];
    locs(locs<4) = [];
    if ~isempty(locs)
        metrics.pk1_time(i) = locs(1);
        metrics.pk1_amp(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time(i) = locs(2);
        metrics.pk2_amp(i) = pks(2);
    end
end
metrics.pk1_time = (metrics.pk1_time-1)/12;
metrics.pk2_time = (metrics.pk2_time-1)/12;


%% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = medfilt1(metrics.time_series,5,[],2);
aux.thresholds = linspace(0,baseline*3,25);
metrics.envelope = zeros(size(metrics.time_series,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*12))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope(i,j)
                    metrics.envelope(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope = metrics.envelope/12;


% Number of frames above a given threshold
metrics.duration = zeros(size(metrics.time_series,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration(:,i) = nansum(smoothed>aux.thresholds(i),2)/12;
end

%% TRIM EVERYBODY to a common length (of "good" sets, our current minimum is about 21 hrs (252 frames)
try
    metrics.time_series = metrics.time_series(:,1:254);
    metrics.integrals = metrics.integrals(:,1:254);
    metrics.derivatives = metrics.derivatives(:,1:252);
catch me
    disp('Note: vectors too short to cap @ 257 frames')
end

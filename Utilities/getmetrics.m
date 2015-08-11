function [metrics,fourier] = getmetrics(id, baseline)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% metrics = getmetrics(id, baseline)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GETMETRICS uses the see_nfkb_native function to filter and preprocess NFkB trajectories,
% then calculates related metrics regarding activation.
% 
% 1) time series (base NFkB dynamics, resampled to 12 frames/hr
% 2) integrated activity
% 3) differentiated activity
% 4) calculated metrics: measuring aspects of oscillation, duration timing ,and amplitude
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Define baseline, if not passed
if nargin <2
    if id <= 270
        baseline = 1;
    else
        baseline = 1.75;
    end
end

% If user passes baseline of 0, reacluclate it from an "off" set (best one is currently #215)
if baseline==0
    [graph] = see_nfkb_native(274); %0.33ng/mL LPS - no activation observed.
    robuststd = @(distr, cutoff) nanstd(distr(distr < (nanmedian(distr)+cutoff*nanstd(distr))));
    baseline = nanmedian(graph.var(:))+2.5*robuststd(graph.var(:),3);
end

% Load/filter/normalize data. Calculate time vector (as function of hrs)
[graph] = see_nfkb_native(id);
t = min(graph.t):1/12:max(graph.t);

% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
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

% 3) differentiated activity - use central fininte difference
smoothed = smoothrows(metrics.time_series,3);
metrics.derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);



% 4) calculated metrics
% GENERAL-PURPOSE: Compute off-time for all cells
metrics.off_times = zeros(size(smoothed,1),1);
window_sz = 43;
tmp = smoothrows(smoothed<(baseline),(window_sz*2)-5);
for j = 1:size(tmp,1)
    thresh_time = find(tmp(j,:)>(1-2/window_sz),1,'first');
    if isempty(thresh_time) 
        thresh_time = find(~isnan(metrics.time_series(j,:)),1,'last');%
    end       
    metrics.off_times(j) = (thresh_time-1)/12;
    if metrics.off_times(j)<0
        metrics.off_times(j) = 0;
    end
end

% METRICS OF OSCILLATION
% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times)*12;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
fourier.fft = zeros(size(metrics.time_series,1),NFFT/2+1);
fourier.freq = Fs/2*linspace(0,1,NFFT/2+1);
fourier.power = zeros(size(fourier.fft));
for i = 1:size(metrics.time_series,1)
    y = metrics.time_series(i,1:depth);
    y((metrics.off_times(i)*12 + 1):end) = nan;
    y(isnan(y)) = [];
    y = y-nanmean(y);
    if ~isempty(y)
        Y = fft(y,NFFT)/length(y);
        fourier.fft(i,:) = abs(Y(1:NFFT/2+1));
        fourier.power(i,:) = abs(Y(1:NFFT/2+1).^2);
    end

end
% Find the point of peak (secondary) power
metrics.peakfreq = nan(size(fourier.power,1),1);
metrics.oscfrac = nan(size(fourier.power,1),1);

for i =1:size(metrics.time_series,1)
    [~,locs] = globalpeaks(fourier.power(i,:),2);
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq(i) = 3600*fourier.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq(i) = 3600*fourier.freq(locs);
    else
        metrics.peakfreq(i) = 3600*fourier.freq(1);
    end
    metrics.oscfrac(i) = nansum(fourier.power(i,(3600*fourier.freq) > 0.5)) /nansum(fourier.power(i,:));
    if isnan(metrics.oscfrac(i))
        metrics.oscfrac(i) = 0;
    end
end




% METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude
metrics.pk1_time = nan(size(smoothed,1),1);
metrics.pk1_amp =  nan(size(smoothed,1),1);
metrics.pk2_time = nan(size(smoothed,1),1);
metrics.pk2_amp =  nan(size(smoothed,1),1);
tmp = smoothed(:,1:96); % restrict to 1st 8 hours

for i = 1:size(tmp,1)
    [amp, time] = globalpeaks(tmp(i,:), 2);
    amp(time<5) = [];
    time(time<5) = [];
    if (~isempty(amp))&&(~isempty(time))
        metrics.pk1_time(i) = min(time);
        metrics.pk1_amp(i) = amp(time==min(time));        
        if length(time)>1
            metrics.pk2_time(i) = max(time);
            metrics.pk2_amp(i) = amp(time==max(time));
        else
            cap_time = max([time+15, depth-1]);
            if ~isempty(nanmax(smoothed(i,cap_time:depth)))
                [metrics.pk2_amp(i), metrics.pk2_time(i)] = nanmax(smoothed(i,cap_time:depth));
            end
            
        end
    
    end

end
metrics.pk1_time = (metrics.pk1_time-1)/12;
metrics.pk2_time = (metrics.pk2_time-1)/12;


% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = medfilt1(metrics.time_series,5,size(metrics.time_series,1),2);

thresholds = baseline/2:(baseline/8):(baseline*2);
metrics.envelope = zeros(size(metrics.time_series,1),length(thresholds));
for j = 1:length(thresholds)
    thresholded = smoothed2>thresholds(j);
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
thresholds = baseline/2:(baseline/8):(baseline*3);
metrics.duration = zeros(size(metrics.time_series,1),length(thresholds));
for i = 1:length(thresholds)
    metrics.duration(:,i) = nansum(smoothed>thresholds(i),2)/12;
end


%% TRIM EVERYBODY to a common length (of "good" sets, our current minimum is about 20 hrs (192 frames)
try
    metrics.time_series = metrics.time_series(:,1:252);
    metrics.integrals = metrics.integrals(:,1:252);
    metrics.derivatives = metrics.derivatives(:,1:250);
catch me
    disp('Note: vectors too short to cap @ 252 frames')
end

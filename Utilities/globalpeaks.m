function [peaks, locs, heights] = globalpeaks(vect, num_peaks)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [peaks, locs] = globalpeaks(vect, num_peaks)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GLOBALPEAKS seeks to find the "dominant" peaks in an input vector - output will be sorted
% from most-to-least dominant.
%
% INPUT:
% vect            input vector
% num_peaks       number of desired "dominant" peaks
%
% OUTPUT:
% peaks          peak values
% locs           peak locations (index)
% heights        peak height (above nearest troughs)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Find all peaks in vector; 1st global peak is maximum value overall
[all_peaks, all_locs] = findpeaks(vect);
all_peaks(((all_locs==1)) | (all_locs==length(vect)) ) = [];
all_locs(((all_locs==1)) | (all_locs==length(vect)) ) = [];

getclosest = @(idx,vect) vect(find(abs(vect-idx)==min(abs(vect-idx)),1,'first'));
peaks = [];
locs  = [];
heights = [];

while length(peaks) < num_peaks
    % Eliminate peaks that have been found already
    tmp_peaks = all_peaks;
    tmp_locs = all_locs;
    tmp_peaks(ismember(tmp_locs,locs)) = [];
    tmp_locs(ismember(tmp_locs,locs)) = [];
    if isempty(tmp_peaks)
      break
    end
    
    % For each candidate, identify nearest peaks - maximize difference btw candidate and two nearest troughs.  
    diffs = zeros(size(tmp_peaks));
    loc_compare = [1 locs length(vect)];
    for i = 1:length(tmp_locs)
        tmp = loc_compare; tmp(tmp>=tmp_locs(i)) = inf;
        trough1 = min(vect(getclosest(tmp_locs(i),tmp):tmp_locs(i)));
        tmp = loc_compare; tmp(tmp<=tmp_locs(i)) = inf;
        trough2 = min(vect(tmp_locs(i):getclosest(tmp_locs(i),tmp)));
        diffs(i) = tmp_peaks(i) - max([trough1, trough2]);
    end
  
    peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
    locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
    heights = [heights, max(diffs)];
end
function [peaks, locs] = globalpeaks(vect, num_peaks)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GLOBALPEAKS seeks to find the "dominant" peaks in an input vector
%
% vect            input vector
% num_peaks       number of desired "dominant" peaks
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

peaks = max(vect);
locs = find(vect==max(vect),1,'first');


while length(peaks) < num_peaks
  [tmp_peaks, tmp_locs] = findpeaks(vect);
  tmp_peaks(ismember(tmp_locs,min(locs):max(locs))) = [];
  tmp_locs(ismember(tmp_locs,min(locs):max(locs))) = [];
  if isempty(tmp_peaks)
      break
  end
  % For each other peak, get the minimum value betweeen it and newest peak.
  diffs = zeros(size(tmp_peaks));
  for i = 1:length(tmp_peaks)
      subset = min([tmp_locs(i),locs(end)]):max([tmp_locs(i),locs(end)]);
      diffs(i) = tmp_peaks(i) - min(vect(subset));
  end
  peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
  locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
end
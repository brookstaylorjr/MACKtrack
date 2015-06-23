function [peaks, locs] = globalpeaks(vect, num_peaks)
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
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Find all peaks in vector; 1st global peak is maximum value overall
[all_peaks, all_locs] = findpeaks(vect);
peaks = max(vect);
locs = find(vect==max(vect),1,'first');

while length(peaks) < num_peaks
  % Eliminate peaks that have been found already
  tmp_peaks = all_peaks;
  tmp_locs = all_locs;
  tmp_peaks(ismember(tmp_locs,locs)) = [];
  tmp_locs(ismember(tmp_locs,locs)) = [];
  if isempty(tmp_peaks)
      break
  end
  % For each candidate peak, identify nearest global peaks - maximize difference btw candidate and the 
  % two surrounding troughs.
  pklist = sort([1 locs length(vect)],'ascend');
  diffs = zeros(size(tmp_peaks));
  for i = 1:length(tmp_peaks)
      begin_ind = pklist(find(pklist<tmp_locs(i),1,'last'));
      end_ind = pklist(find(pklist>tmp_locs(i),1,'first'));
      diffs(i) = tmp_peaks(i) - (min(vect(begin_ind:tmp_locs(i))) + min(vect(tmp_locs(i):end_ind)))/2;
  end
  peaks = [peaks, tmp_peaks(find(diffs==max(diffs),1,'first'))];
  locs = [locs, tmp_locs(find(diffs==max(diffs),1,'first'))];
end
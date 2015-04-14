function [thresh, noiseCount, val] =noisethresh(image_in,dropPixels,searchRange,noiseSize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [thresh, noiseCount, val] =noisethresh(image_in,dropPixels,searchRange,noiseSize)
%
% NOISETHRESH computes an image threshold based on the assumption that optimal placement
% is at the LOWEST value that still minimizes extra background noise.
%
% image_in      input image (intensity or edge magnitudes)
% dropPixels    background area of images; normalizes over cell densities
% searchRange   min and max values that set threshold search range 
% noiseSize     maximum size of "noisy" pixel groups
%
% thresh        output threshold
% noiseCount    number of noise-sized object pixels at each value tested
% val           edge values tested
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Take input image and scan over values, constructing noiseCount vector
val = min(searchRange): (max(searchRange)-min(searchRange))/100:max(searchRange);
image_in(dropPixels) = 0;
noiseCount = zeros(size(val));

for i = 1:length(val)
    mask = (abs(image_in)>val(i));
    mask2 = bwareaopen(mask,noiseSize);
    sum1 = sum(mask(:)-mask2(:));
    noiseCount(i) = sum1;
end

try
    % 2. Find difference for function
    noiseDiff = noiseCount(3:end) - noiseCount(1:end-2); % Get derivative (central difference) of noiseCount
    
    % 3. Smooth noiseDiff and calculate "elbow" slope
    indNeg = find(noiseDiff==min(noiseDiff),1,'first'); % Most negative slope
    elbowSlope = (noiseDiff(indNeg))*0.33; % "Elbow" slope: ~1/3 of most negative slope.
    
    % 4. Find first time slope drops beneath elbowSlope for 2 consecutive pts
    flag = 0;
    ind = indNeg;
    noiseDiff_1 = noiseDiff;
    noiseDiff_1(1:ind) = nan;
    while (~flag) && (sum(isnan(noiseDiff_1))~=sum(isnan(noiseDiff)))
        noiseDiff = noiseDiff_1;
        ind = find(noiseDiff > elbowSlope, 1,'first');
        if (sum(noiseDiff(ind:ind+4)>elbowSlope)>3)
            flag = 1;
        end
        noiseDiff_1(1:ind) = nan;
    end
    thresh = val(ind);
    
catch exception
    disp('Edge thresholding failed; threshold set to midpoint of search range.')
    thresh = mean(searchRange);
    noiseCount = [];
end

if isempty(thresh)
    thresh = mean(searchRange);
    noiseCount = [];
end
function [matches] = matchclosest(candidates, seeds,imageSize)
% MATCHCLOSEST Match a single point from a group of 'candidates' to its nearest 'seed' location.
%
% candidates    all candidate locations:  can be binary mask or a 2-row array of indicies [rows; cols]
% seeds         all seed locations:  can be binary mask or a 2-row array of indicies [rows; cols]
% imageSize     optional input, changes form of output matches
%
% matches       one match per seed is output (in row/col array), unless imageSize is defined- in that case, matches is a binary
% array         showing chosen candidates.

% Determine format of data: if the data is a 2xn array, keep as is. Otherwise, assume it's a binary mask and format correctly.
candidates(isnan(candidates)) = 0;
seeds(isnan(seeds)) = 0;

if size(candidates,1) > 2
    if max(candidates(:))>1
        error('Input data not formatted correctly')
    else
    [row, col] = find(candidates==1);
    candidates = [row,col]';
    end
end

if size(seeds,1) > 2
    if max(seeds(:))>1
        error('Input data not formatted correctly')
    else
    [row, col] = find(seeds==1);
    seeds = [row,col]';
    end
end

% Make structure to hold the output
matches = zeros(size(seeds));
for i = 1:size(seeds,2)
    dist = ( (candidates(1,:)-seeds(1,i)).^2 + (candidates(2,:)-seeds(2,i)).^2 );
    matches(:,i) = candidates(:,find(dist == min(dist),1,'first'));
end

if nargin>2
    matchesInd = sub2ind(imageSize,matches(1,:),matches(2,:));
    matches = false(imageSize);
    matches(matchesInd) = 1;
end
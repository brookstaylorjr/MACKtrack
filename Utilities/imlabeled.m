function [] = imlabeled(imageDisp, labelMat, dispString, imgZoom, textColor)
% IMLABELED Display selected object characteristic at the centroid of each object in a label matrix.
% 
% INPUTS
% -imageDisp: the base image to be displayed
% -labelMat: a label structure/matrix, used to calculate centroids for display (can also be a binary mask)
% -dispString: the text to be displayed: can be a single character, 1xn string cell array, or 'Area' (calc. from regionprops)
% -imgZoom: size of image- 0.85 is good for 1024x1248 images
% -textColor:  triplet denoting color- defaults to white ([1 1 1])


% Get text color
if nargin<5
    textColor = [1 1 1];
end
% Parse the dispString input- get centroids (and if applicable, areas)
if ischar(dispString)
    if strcmp(dispString,'Area') 
        s = regionprops(labelMat,'Centroid',dispString);
        dispString = cell(size(s));
        for i = 1:length(s)
            dispString{i} = num2str(s(i).Area);
        end
    end
else
    s = regionprops(labelMat,'Centroid');
end

% Parse centroids
obj = unique(labelMat);
obj(obj ==0) = [];

markerX = zeros(size(obj));
markerY = zeros(size(obj));
for i = 1:length(obj)
    markerX(i) = s(obj(i)).Centroid(1);
    markerY(i) = s(obj(i)).Centroid(2);
end
        
        
figureLoc = [10,45,round(size(imageDisp,2)*imgZoom),round(size(imageDisp,1)*imgZoom)];
set(gcf,'Position',figureLoc);
imageDisp = double(imageDisp);
imshow(imageDisp/max(imageDisp(:)),'Border','tight','InitialMagnification','fit')

text(markerX-5,markerY+5,dispString,'FontSize',7,'Color',textColor)

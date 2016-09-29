function [] = saveFig(inputImg,CellLabel,NuclearLabel, X, bit_depth, saveName, alpha, outputRes, saturationVal)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SAVEFIG Save high-quality display figure showing pseudocolor nuclei on top of outlined cells
%
% [] = saveFig(inputImg,CellLabel,NuclearLabel, X, bit_depth, saveName, alpha, outputRes, saturationVal)
%
%
% inputImg       base image (DIC/phase) to be saved
% CellLabel      cell label matrix
% NuclearLabel   nuclear label matrix
% X              background correction matrix (pass empty vector, if no correction)
% bit_depth      bit depth of base image
% saveName       output file name
% alpha          transparency for overlay - 0.3 works for transmitted light images (default), 
%                    and < 0.2 works better for fluorescence
% outputRes      output image resolution (default = [768 1024])
% saturationVal  range of output image (img is rescaled so background = 0, with stdev of 1) - default = [-2 4]
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - Setup: define output resolution, rotate images if necessary - - - -
if nargin < 9
	saturationVal = [-2 4];
end
if nargin < 8
    max_dim = max(size(inputImg));
    outputRes = round([size(inputImg,1)/max_dim*1024 size(inputImg,2)/max_dim*1024]);
end
if nargin < 7
    alpha = 0.3;
end


% Check if image is taller than it is wide- rotate it 90 deg if that's the case
if size(inputImg,1)>size(inputImg,2)
   inputImg = imrotate(inputImg,90);
   CellLabel = imrotate(CellLabel,90);
   NuclearLabel = imrotate(NuclearLabel,90);
   outputRes = [outputRes(2),outputRes(1)];
end

% - - - - Make base phase image: drop minimum to zero, balance mode at 90, and convert to 8-bit. - - - -

% Background correct image (if necessary)
if ~isempty(X)
    % Calculate background function
    warning off MATLAB:nearlySingularMatrix
    pStar = (X'*X)\(X')*double(inputImg(:));
    % Apply correction
    inputImg = reshape((double(inputImg(:) - X*pStar)),size(inputImg));
    inputImg = inputImg-min(inputImg(:)); % Set minimum to zero
end

% Balance background mode of image to 0, then map to 0-255
baseImg = modebalance(inputImg,1,bit_depth);
x = saturationVal;
baseImg(baseImg<x(1)) = x(1);
baseImg(baseImg>x(2)) = x(2);
baseImg = (baseImg - x(1))/(x(2)-x(1))*255;
baseImg = uint8(round(baseImg));

% - - - - Make overlay image: define colormap for cells/borders, convert label matrix  - - - -
borders= ((imdilate(CellLabel,ones(3))-CellLabel)>0);
map = [0 81 84;
    19 41 92;
    76 19 140;
    84 0 7;
    28 113 42;
    30 65 161;
    60 61 64;
    48 88 140;
   240 157 0]/255;
numColors = size(map,1);
cellLabel = (CellLabel>0).*(mod(double(CellLabel),numColors-1)+1);
cellLabel(borders) = numColors;
labelRGB = label2rgb(cellLabel,map,'w');
% define transparency
transparency = alpha.*(double(baseImg)/255);
transparency(borders) = transparency(borders)*1.9;
transparency(transparency>1) = 1;
transparency(cellLabel==0) =  0;

% composite baseImg together with labelRGB
composite = cat(3, uint8(double(baseImg).*(1-transparency) + double(labelRGB(:,:,1)).*transparency),...
            uint8(double(baseImg).*(1-transparency) + double(labelRGB(:,:,2)).*transparency),...
            uint8(double(baseImg).*(1-transparency) + double(labelRGB(:,:,3)).*transparency));

% resize image, get scale factor so centroids aren't off
composite = imresize(composite,outputRes);

% - - - - Add text: Get nuclear centroids and label with cell ID number, then add annotations - - - -
s = regionprops(NuclearLabel,'Centroid');
obj = unique(NuclearLabel(NuclearLabel>0));
try % (add fallback for older versions of MATLAB)
    locs = cell2mat(struct2cell(s(obj)));
    locs = [(locs(1:2:end))',(locs(2:2:end))'];
    composite = insertText(composite,locs,obj,'FontSize',12,'TextColor',[255,255,255],...
        'Font','Arial','BoxOpacity',0,'AnchorPoint','center');
catch me
    scale_factor = size(composite,1)/size(NuclearLabel,1);
    nums =  cellstr(num2str(obj(:)))';
    txtInserter = vision.TextInserter(nums,'LocationSource','Input port','FontSize', 12,'Color',[255,255,255]);

    for i = 1:length(obj)
        loc = uint16(s(obj(i)).Centroid*scale_factor - [7 7]);
        composite = step(txtInserter,composite,uint16(i),loc);
    end
end
% Write image to file 
 imwrite(composite,saveName,'jpg','Quality',90)
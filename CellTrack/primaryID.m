function [output, diagnos] = primaryID(nucOrig,p,~)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PRIMARYID:     Identify nuclei in fluorescent image (Hoechst stain) for NF-kB tracing
%
% nucOrig        original fluorescent nuclear image
% p              parameters struture
% X              background correction matrix
%
% label_final    output label matrix with nuclei
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
min_area = round(pi*(p.MinNucleusRadius^2));
max_area = round(pi*(p.MaxNucleusRadius^2));
cutoff.Area = [min_area, max_area];
cutoff.Eccentricity = p.Eccentricity;
cutoff.Compactness = p.Compactness;
diagnos.nucleus1 = nucOrig;


%- - - - - - - - - - - - - - - - - - - MASK1 - - - - - - - - - - - - - - - - - - - - - - -
% Initial Otsu threshold on image
diagnos.thresh1  = otsuthresh(diagnos.nucleus1,false(size(diagnos.nucleus1)),'none'); 
% Log compression is strongly affected by object content of image
diagnos.mask1a = diagnos.nucleus1>diagnos.thresh1;
% Cleanup
diagnos.mask1b =  bwareaopen(diagnos.mask1a,p.NoiseSize,4); % Despeckle
diagnos.mask1b = ~bwareaopen(~diagnos.mask1b,256,4); % Hole fill
diagnos.mask1b = imopen(diagnos.mask1b,strel('diamond',3));
diagnos.mask1b = bwareaopen(diagnos.mask1b,min_area);

% Break up shapes by distance transform
dist_smoothed = bwdist(~diagnos.mask1b);
dist_smoothed = imfilter(dist_smoothed,gauss2D(p.MinNucleusRadius/4),'replicate');
dist_smoothed = medfilt2(dist_smoothed,[5 5]);
diagnos.label1a = imdilate(watershedalt(dist_smoothed,diagnos.mask1b,4),ones(3));

% Try to bridge back watershed-split nuclei, then remove all small objects
diagnos.label1 = bridgenuclei(diagnos.label1a,cutoff);

%- - - - - - - - - - - - - - - - - - - MASK2 - - - - - - - - - - - - - - - - - - - - - - -
diagnos.nucleus2 = imfilter(diagnos.nucleus1,gauss2D(p.MinNucleusRadius/4),'replicate'); % Gaussian filtered
diagnos.thresh2  = otsuthresh(diagnos.nucleus2,false(size(diagnos.nucleus2)),'log');
diagnos.mask2a = diagnos.nucleus2>diagnos.thresh2;
diagnos.mask2b = ~bwareaopen(~diagnos.mask2a,256,4); % Hole fill
diagnos.mask2b = imopen(diagnos.mask2b,diskstrel(round(p.MinNucleusRadius/5)));
diagnos.mask2b = imclose(diagnos.mask2b,diskstrel(round(p.MinNucleusRadius/4))); % Do larger close here

% Break up mask by distance matrix
dist_smoothed = bwdist(~diagnos.mask2b);
dist_smoothed = imfilter(dist_smoothed,gauss2D(p.MinNucleusRadius/4),'conv');
dist_smoothed = medfilt2(dist_smoothed,[p.MedianFilterSize,p.MedianFilterSize]);
diagnos.label2a = watershedalt(dist_smoothed,diagnos.mask2b,4);
markers = bwmorph(diagnos.label1>0,'shrink',15);
diagnos.label2b = imdilate(removemarked(diagnos.label2a,markers,'remove'),ones(3));

% Try to bridge back watershed-split nuclei, then remove all small objects
diagnos.label2 = bridgenuclei(diagnos.label2b,cutoff);

% Combine label matricies
diagnos.label2(diagnos.label1>0) = 0; % Double check and make sure there's no overlap
unique1 = unique(diagnos.label1);
unique1(unique1==0) = [];
for i = 1:length(unique1)
    diagnos.label1(diagnos.label1==unique1(i)) = i;
end

unique2 = unique(diagnos.label2);
unique2(unique2==0) = [];
start_val = max(diagnos.label1(:));
for i = 1:length(unique2)
    diagnos.label2(diagnos.label2==unique2(i)) = i;
end
diagnos.label2(diagnos.label2>0) = diagnos.label2(diagnos.label2>0) + start_val;
label_final = diagnos.label1+diagnos.label2;
label_final = imclose(label_final,(diskstrel(round(p.MinNucleusRadius/3))));
label_final = imopen(label_final,(diskstrel(round(p.MinNucleusRadius/2))));
output.nuclei = label_final;

% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);
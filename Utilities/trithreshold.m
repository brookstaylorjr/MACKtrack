function [thresh, all_thresh] = trithreshold(img,bins,mode_level)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [thresh, all_thresh] = trithreshold(img,bins,mode_level)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% TRITHRESHOLD performs 3 different thresholding methods on an image:
% 1) Otsu's method (maximizes inter-population variance)
% 2) Tsai's method (maximizes a change in curvature on the histogram graph)
% 3) Gaussian mixture modeling
%
% INPUTS:
% img           input image (should have distinct background and foreground regions)
% bins          bins to divide image into (used for MoG modeling)
% mode_level   mode of image (1st guess is that this represents background in a GMM)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

v = img(:);
v(v==max(v(:))) = [];
v(v==min(v(:))) = [];

if nargin<2
        bins = linspace(prctile(v,0.5),prctile(v,99.5),512);
end

if nargin<3
    N = histcounts(v,bins);
    N(1) = 0; N(end) = 0;
    mode_level = bins(find(N==max(N),1,'first'));
end



all_thresh = zeros(3,1);
% 1) Otsu threshold
all_thresh(1) = quickthresh(img,false(size(img)),'none');

% 2) Tsai threshold (based on histogram shape)
try
    all_thresh(2) = tsaithresh(v);
catch me
    all_thresh(2) = nan;
end

% 3) Mixture model of 2 Gaussian distributions
start_2p = struct;
sigma1 = std(v(v<mode_level))*5/3; % Get std from bottom half of left distribution
start_2p.mu = [mode_level; mode_level+3*sqrt(sigma1)];
start_2p.Sigma = cat(3,sigma1, (sigma1*3));
p = min([0.97 numel(v(v<mode_level))*2/numel(v)]);
start_2p.PComponents = [p 1-p];
% Set maxIter for GM modeling (should be small for speed, <50)
gmopt = statset('MaxIter',50);
% Model distribution - subsample image a little bit (100K points)
warning('off','stats:gmdistribution:FailedToConverge')
dist_2p = fitgmdist(v(1:4:end),2,'Start',start_2p,'CovType','diagonal','Options',gmopt);
warning('on','stats:gmdistribution:FailedToConverge')
tmp = dist_2p.cluster(bins');
tmp_switch = [1+find(abs(diff(tmp))>0)'];
if length(tmp_switch)>1
    len1 = [tmp_switch(1)-1 length(tmp)-tmp_switch(2)];
    tmp_switch = tmp_switch(find(len1==max(len1),1,'first'));
end
try
all_thresh(3) = bins(tmp_switch);
catch me
    all_thresh(3) = nan;
    warning(['GMM thresholding failed (output was [',num2str(bins(tmp_switch)),'])- omitting.'])
end
thresh = nanmedian(all_thresh);
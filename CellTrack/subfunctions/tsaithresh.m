function [thresh,K,H,bins,testFn] = tsaithresh(image1,dropPixels,numBins)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%  TSAITHRESH Threshold image using method by Tsai (1995)
%
% Inputs
% - image1: image to be thresholded
% - dropPixels: binary mask showing pixels to be pulled out of the image
%
% Outputs
% - thresh: calculated image threshold, above mode of image
% - K : curvature of smoothed histogram
% - H: smoothed histogram function
% - bins: histogram bins from imagefrom
%
% NOTE: we are dealing with phase images, which are generally unimodal- the
% multimodal part of this algorithm is currently omitted
%
% Tsai, D.M. "A fast thresholding selection procedure 
% for multimodal and unimodal histograms"
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<3
    numBins = 4096;  % Default; good for 12/14-bit images
end

% Get image histogram
image1(dropPixels) = [ ];
image1 = image1(:);
image1(image1==max(image1)) = [];

% Default parameters, should not need to be changed.
fsize = 3; % Gaussian filter size
R =2; % Radius of curvature, Tsai says <= 3
lowCutoff = 0.02; % sets cutoff for very low peaks

[histFunc, bins] = hist(image1,numBins);


% Make smoothing filter
gaussian = [0.2261 0.5478 0.2261];
gauss = gaussian;
if fsize > 1
    for i = 1:fsize-1
        gauss = conv(gauss,gaussian);
    end
end

% Smooth function, find peaks- Drop all peaks that are less than 2% of max.
 H = ifft(fft(histFunc).*fft(gaussian,length(histFunc)).*conj(fft(gaussian,length(histFunc))));
 [peaks locs] = findpeaks(H, 'minpeakheight',max(H)*lowCutoff);
 index = 0;


% Continue to smooth til there is only one large peak
 while (length(peaks) > 1) && (index<10)
    H = ifft(fft(H).*fft(gaussian,length(H)).*conj(fft(gaussian,length(H))));
 [peaks locs] = findpeaks(H, 'minpeakheight',max(H)*lowCutoff);
 index = index+1;
 end


% ____Unimodal Distribution___________________________________________________
%(maximize rate of change of curvature)
    
% Calculate psi  and K at each point for the histogram function
tCalc = 1+ R*2 : length(H) - R*2; % indices that we're going to check over
psi = zeros(size(H));
K = zeros(size(H));
for  t = tCalc
    for j = 1:R
        psi(t) = psi(t) + (1/R).*(H(t+j) - H(t-j))/(2*j);
    end
end
for  t = tCalc
    for j = 1:R
        K(t) = K(t) + (1/R).*abs(psi(t+j) - psi(t-j));
    end
end


 % Start looking for a maxima after peak H/ peak K value (whichever is farther)
[Kpeaks Klocs]= findpeaks(K);
ind1 = max([find(Kpeaks==max(Kpeaks)), find(Klocs<=find(H==max(H)))]);
Kpeaks_drop = Kpeaks;
Kpeaks_drop(1:ind1) = [ ];
ind2 = find(Kpeaks == Kpeaks_drop(1));
testFn = K(find(K==Kpeaks(ind2-1)):find(K==Kpeaks(ind2)));

% Need to watch out for false positive on the way down from previous peak.
% check: if VALLEY between highest peak and next highest peak is also high,
% find the next highest.
counter = 0;
break_flag = 0;
while (min(testFn) > (Kpeaks_drop(1)/2)) && (counter<3)
    try
        Kpeaks_drop(1) = [];
        ind2 = ind2+1;
        testFn = K(find(K==Kpeaks(ind2-1)):find(K==Kpeaks(ind2)));
        counter = counter+1;
    catch me
        break_flag = 1;
        break
    end
end

if ~break_flag
    thresh = bins(K==Kpeaks_drop(1));
else
    thresh = otsuthresh(image1,dropPixels,'none');
end


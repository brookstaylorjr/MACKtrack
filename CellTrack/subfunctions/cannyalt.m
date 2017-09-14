function edgemask = cannyalt(image0, gauss_size)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% edgemask = cannyalt(image0, gauss_size)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CANNYALT performs Canny edge-finding with alternate (Tsai) thresholding on data - this 
% threshold tends to perform more robustly than MATLAB default.
%
% image0         DIC image
% gauss_size     std dev of gaussian filter (diameter is created automatically)
%
% edgemask     output mask with thinned cell edges
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Gaussian-filter image to find edges at lower resolutions 
imfilt1 = imfilter(image0, gauss2D(gauss_size),'replicate');

% Calculate edges, then magnitude and direction
horizontalEdge = imfilter(imfilt1,fspecial('sobel') /8,'replicate');
verticalEdge = imfilter(imfilt1,fspecial('sobel')'/8,'replicate');
edge_mag = sqrt(horizontalEdge.^2 + verticalEdge.^2);
edge_dir = atan(horizontalEdge./verticalEdge);

% Round edge direction to nearest 45 degrees
edge_0 = (edge_dir >=  -pi/8) & (edge_dir < pi/8);
edge_45 = (edge_dir >= (pi/8)) & (edge_dir < (3*pi/8));
edge_90 = abs(edge_dir) >= (3*pi/8);
edge_135 = (edge_dir < (-pi/8)) & (edge_dir > -(3*pi/8));

% Thin edges by comparing each pixel to its immediate neighbor
edge_thin = (edge_0(2:end-1,2:end-1) & ((edge_mag(2:end-1,2:end-1) > edge_mag(2:end-1,1:end-2)) & (edge_mag(2:end-1,2:end-1) > edge_mag(2:end-1,3:end)))) |...
(edge_45(2:end-1,2:end-1) & ((edge_mag(2:end-1,2:end-1) > edge_mag(1:end-2,1:end-2)) & (edge_mag(2:end-1,2:end-1) > edge_mag(3:end,3:end)))) |...
(edge_135(2:end-1,2:end-1) & ((edge_mag(2:end-1,2:end-1) > edge_mag(3:end,1:end-2)) & (edge_mag(2:end-1,2:end-1) > edge_mag(1:end-2,3:end)))) |...
(edge_90(2:end-1,2:end-1) & ((edge_mag(2:end-1,2:end-1) > edge_mag(1:end-2,2:end-1)) & (edge_mag(2:end-1,2:end-1) > edge_mag(3:end,2:end-1))));
% Pad out thinned edge image to full size again
edge_thin = [false(size(edge_thin,1),1),edge_thin,false(size(edge_thin,1),1)];
edge_thin = [false(1,size(edge_thin,2));edge_thin;false(1,size(edge_thin,2))];

% Low threshold: Tsai threshold
[t1,~,H,bins] = tsaithresh(edge_mag,~edge_thin);
edgemask_low = (edge_mag>t1) & edge_thin;

% High threshold: look at local slope of smoothed hist function (H) ->  
% find point where slope decreases to 10% of value @ Tsai threshold
H_slope = diff(H);
x = bins(1:end-1);
H_slope(x<t1) = [];
x(x<t1) = [];
t2 = x(find(abs(H_slope) < 0.10*(abs(H_slope(1))),1,'first'));

% Fallback: if Tsai threshold did a bad job (e.g. function was very "bumpy"), I need a different fallback -> 
% find point of max negative slope, then keep going til you hit 10% of that value
if isempty(t2)
    H_slope(1:find(H_slope==min(H_slope),1,'first')) = inf;
    
    t2 = x(find(abs(H_slope) < 0.10*max(abs(H_slope)),1,'first'));
end
edgemask_high = (edge_mag>t2) & edge_thin;

% Hysteresis thresholding: use removemarked to only keep strong-marked edgea
edgemask = labelmatrix(removemarked(bwconncomp(edgemask_low,8),edgemask_high,'keep'))>0;


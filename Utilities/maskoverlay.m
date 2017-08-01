function RGB = maskoverlay(img_in, mask, overlay_color, a)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% RGB = maskoverlay(img_in, mask, overlay_color, a)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MASKOVERLAY creates an RGB composite image of (a) a base image (can be RGB or BW) with
% (b) an overlaid mask, with specified transparency, a
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<4
    a = 0.5;
end

if ismatrix(img_in)
    R = img_in; G= img_in; B = img_in;
else
    R = img_in(:,:,1); G = img_in(:,:,2); B = img_in(:,:,3);
end

R(mask) = img_in(mask)*(1-a) + overlay_color(1)*a;
G(mask) = img_in(mask)*(1-a) + overlay_color(2)*a;
B(mask) = img_in(mask)*(1-a) + overlay_color(3)*a;
RGB = cat(3,R,G,B);
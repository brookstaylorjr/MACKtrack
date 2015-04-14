function SE = diskstrel(radius)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DISKSTREL constructs a disk-shaped structuring element (via 'strel'), but chooses the
% number of periodic elements used in the approximation of the circle.
% 
% radius    input radius of the sturcturing element
%
% SE        output structuring element
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


if (radius<= 8) 
    n = 4;
elseif (radius<=11)
    n = 6;
elseif (radius<=16)
    n = 8;
else
    n = 0;
end
SE = strel('disk',radius,n);
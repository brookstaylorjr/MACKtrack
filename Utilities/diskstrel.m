function SE = diskstrel(radius)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SE = diskstrel(radius)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DISKSTREL automatically generates a nice, round disk-shaped structuring element (via 'strel'),
% of specified radius.
%
% INPUT:
% radius    input radius of the sturcturing element
%
% SE        output structuring element
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

radius = round(radius);
if radius<=0
	error('diskstrel radius must be positive')
end

% Small radii (<=6)
switch radius
    case 1
        SE = strel('disk',radius,8);
    case 2
        SE = strel('disk',radius,8);
    case 3
        SE = strel('octagon',radius);
    case 4
        SE = strel('disk',radius,0);
    case 5
        SE = strel('disk',radius,4);
    case 6
        SE = strel('octagon',radius);

end

% Larger radii
if (radius>6) && (radius<=8)
    SE = strel('disk',radius,4);
elseif (radius>8) && (radius <=11)
    SE = strel('disk',radius,6);
elseif (radius>11) && (radius<=16)
    SE = strel('disk',radius,8);
elseif (radius>16)
    SE = strel('disk',radius,0);
end


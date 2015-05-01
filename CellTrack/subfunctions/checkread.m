function img = checkread(name, bitdepth, double_flag,debug)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [img ] = CHECKREAD(name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CHECKREAD reads in an image, ensuring that it isn't max-padded
%
% INPUT:
% name          name of image (+ file path)
% bitdepth      bit depth of image
% double_flag   if 1, outputs double of image (as opposed to read precision) (default=1)
% debug         if 1, show initial/final image size
%
% OUTPUT:
% img       output image (double)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<3
    debug = 0;
end
if nargin<3
    double_flag = 1;
end
if nargin<2
    bitdepth = 16;
end
max_val = (2^bitdepth)-1;
%min_val = 0;

img = imread(name);
old_nrow = size(img,1);
old_ncol = size(img,2);

% Identify if there are any saturated columns/rows padded at edges of image - drop, if present
ind_r = sum(img==max_val,2)>=(old_ncol*0.99);
if sum(ind_r)>0
    if (find(ind_r,1,'first')==1) || (find(ind_r,1,'last')==old_nrow)
        img(ind_r,:) = [];
    end
end

ind_c = sum(img==max_val)>(old_nrow*0.99);
if sum(ind_c)>0
    if (find(ind_c,1,'first')==1) || (find(ind_c,1,'last')==old_ncol)
        img(:,ind_c) = [];
    end
end
%img(sum(img==min_val,2)>=(old_nrow*0.99),:) = [];
%img(:,sum(img==min_val)>(old_ncol*0.99)) = [];


new_nrow = size(img,1);
new_ncol = size(img,2);

if debug
    if (old_nrow~=new_nrow)||(old_ncol~=new_ncol)
        disp(sprintf(['Note: max-saturated rows/cols dropped from "', name,'"\n',...
                'Old dimensions: [',num2str(old_nrow),', ',num2str(old_ncol),']. ',...
                'New dimensions: [',num2str(new_nrow),', ',num2str(new_ncol),']']));
    end
end

if double_flag
    img = double(img);
end
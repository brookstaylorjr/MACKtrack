function X = backgroundcalculate(imgsize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BACKGROUNDCALCULATE Used in applying flatfield correction (method from Li et al (2008))
%
% imgsize      size ([r,c]) of input image to be corrected
%
% X            output correction matrix
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Li K, Kanade T. Nonnegative mixed-norm preconditioning for microscopy image
% segmentation. Information processing in medical imaging. 2009;21(DIC): 362-73. 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

K = 3;
X = zeros(imgsize(1)*imgsize(2),(K+1)*(K+2)/2,'single');
for i = 1:imgsize(2)
    for j = 1:imgsize(1)
       X(j+(i-1)*imgsize(1),:) = single([1 i j i.^2 i.*j j.^2 i.^3 j.*i.^2 i.*j.^2 j.^3 ]);

    end
end


function X = backgroundcalculate(imgsize,K)
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

if nargin<2
    K=3;
end

X = zeros(imgsize(1)*imgsize(2),(K+1)*(K+2)/2,'single');


switch K
    case 2
        for i = 1:imgsize(2)
            for j = 1:imgsize(1)
               X(j+(i-1)*imgsize(1),:) = single([1 i j i.^2 i.*j j.^2 ]);

            end
        end
    case 3
        for i = 1:imgsize(2)
            for j = 1:imgsize(1)
               X(j+(i-1)*imgsize(1),:) = single([1 i j i.^2 i.*j j.^2 i.^3 j.*i.^2 i.*j.^2 j.^3 ]);

            end
        end
    case 4
        for i = 1:imgsize(2)
            for j = 1:imgsize(1)
               X(j+(i-1)*imgsize(1),:) = single([1 i j i.^2 i.*j j.^2 i.^3 j.*i.^2 i.*j.^2 j.^3 ...
                   i.^4 i.^3*j i.^2*j.^2 i*j.^3 j.^4]);

            end
        end
    otherwise
        error('Flatfield bias can only be fit with quadratic (K=2), cubic (K=3) or quartic (K=4) function.')

end

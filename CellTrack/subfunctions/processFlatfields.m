function [flatfields_out] = processFlatfields(flatfields_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [flatfields_in] = processFlatfields(flatfields_in, image_size)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PROCESSFLATFIELDS takes in a cell matrix of flatfields, and:
% - fits them with a 2-D cubic function (best fits existing measured flatfields)
% - subtracts a darkfield image
% - appends the darkfield image to the existing cell matrix
%
% NOTE: existing darkfield images are stored in CellTrack/subfunctions - need to reflect YOUR
% scope's setup (at the moment, corresponding darkfield can be selected by image size).
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


if isempty(flatfields_in)
    flatfields_out = {};
    return;
end

sz = size(flatfields_in{1});


% a) Grab the correct darkfield image (done based on image's size...)
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end));
BG = imread([home_folder,filesep,'darkfield_',num2str(sz(1)),'x',num2str(sz(2)),'.tif']);



% b) Calculate flatfield images - replace them in parameters
X = backgroundcalculate(sz,3);
warning off MATLAB:nearlySingularMatrix
for i = 1:length(flatfields_in)
    corr_img = double(flatfields_in{i}) - double(BG);
    if min(corr_img(:)) < 0
        error('Flatfield normalization yielded negative values. Please reload these images.')
    end
    
    pStar = (X'*X)\(X')*corr_img(:);
    % Apply correction
    corr_img = reshape(X*pStar,size(corr_img));
    flatfields_in{i} = corr_img;
end
warning on MATLAB:nearlySingularMatrix

% c) Append darkfield image
flatfields_out = [flatfields_in,{BG}];
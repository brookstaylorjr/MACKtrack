function [output_img, varargout] = modebalance(input_img, num_modes, bit_depth, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% modebalance takes any input image and attempts to normalize the background distribution,
% i.e. the first "major" peak. Used in making output images and fluorescence measurements.
%
% [output_img, varargout] = modebalance(input_img, num_modes, bit_depth, varargin)
% input_img      raw input image (can be any imaging modality)
% num_modes      number of subpopulations (should be one or two)
% bit_depth      bit depth of image (used in assigning number of bins)
%
% output_img     adjusted image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Validate number of subpopulations (needs to be 1 or 2)
if ~ismember(num_modes,[1 2]);
    error('Error in MODEBALANCE: only uni- or bi-modal populations can be modeled')
end

% Prep image: drop saturated pixels, convert to vector
imgvect = input_img(:);
imgvect(imgvect==min(imgvect)) = [];
imgvect(imgvect==max(imgvect)) = [];

% Generate histogram from imgvect
x = 0:4:(2^bit_depth)-1; 
n = hist(imgvect,x)/numel(imgvect);
% Get mode of distribution, and calculate percentiles on left side
mode_loc = find(n==max(n),1,'first');
pct1 = cumsum(n(mode_loc:-1:1)) / sum(n(1:mode_loc));
% Determine 1st guess at background mean and variance
mu1 = x(mode_loc);
if max(pct1>0.68)
    sigma1 = (x(mode_loc) - x(mode_loc - find(pct1>0.68,1,'first')))^2;
else
    sigma1 = x(mode_loc).^2;
end


% Turn convergence warning off
warning('off','stats:gmdistribution:FailedToConverge')

switch num_modes
    case 1 % Unimodal case
        % Assign starting GMM structure
        S.mu = mu1;
        S.Sigma = sigma1;
        S.PComponents = 1;
        % Set maxIter for GM modeling (should be small for speed, <50)
        gmopt = statset('MaxIter',5);
        % Model distribution
        obj = gmdistribution.fit(imgvect,1,'Start',S,'SharedCov',true, 'Options',gmopt);  
    case 2 % Bimodal case
        % Assign starting GMM structure
        S.mu = [mu1; mu1+2*sqrt(sigma1)];
        S.Sigma = cat(3,sigma1, sigma1*4);
        S.PComponents = [0.75 0.25]; % This may not be the best thing to do, here...
        % Set maxIter for GM modeling (should be small for speed, <50)
        gmopt = statset('MaxIter',12);
        % Model distribution
        obj = gmdistribution.fit(imgvect,2,'Start',S,'CovType','diagonal','Options',gmopt);
end

% Turn convergence warning back on
warning('on','stats:gmdistribution:FailedToConverge')


% Get the minimum mean and its standard deviation (background) and normalize that to [0,1]
mu_ind = find(obj.mu == min(obj.mu));
background_mu = obj.mu(mu_ind);
if size(obj.Sigma,1) > 1
    background_sigma = sqrt(obj.Sigma(mu_ind,mu_ind));
else
    background_sigma = sqrt(obj.Sigma(mu_ind));
end

% 3 possible modes: 
% - 'display'    (normalize background to 0, w/ standard deviation of 1) - defualt
% - 'measure'    return background distribution only
% - 'correct'    normalize background to match an input distribution
if isempty(varargin)
    varargin{1} = 'display';
end

if strcmp(varargin{1},'measure')
    varargout{1} = [background_mu(1); background_sigma(1)];
    output_img = input_img; % Don't modify image
elseif strcmp(varargin{1},'correct')
    if length(varargin)<2
        error('Error in MODEBALANCE: corrected images need a distribution to match of form [mu; sigma]')
    end
    match_dist = varargin{2};
    output_img = (input_img-background_mu(1))/background_sigma(1)*match_dist(2); % Match standard deviation
    output_img = output_img + match_dist(1); % Match mean
    output_img(output_img<0) = 0;
    output_img(output_img > ((2^bit_depth)-1)) = (2^bit_depth)-1;
else % display
    output_img = (input_img - background_mu(1))/(background_sigma(1));
end



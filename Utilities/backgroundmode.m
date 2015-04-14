function [background_mu] = backgroundmode(image_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BACKGROUNDMODE calculates background mode of a particular image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Get background mu of both NFkB and nuclear images
% Set options
numIter = 5;
gmopt = statset('MaxIter',numIter);
% Drop saturated pixels
imgvect = image_in(:);
imgvect(imgvect==min(imgvect)) = [];
imgvect(imgvect==max(imgvect)) = [];
% Turn (expected) warning off
warning('off','stats:gmdistribution:FailedToConverge')
% We gonna try modeling this thang as a bimodal population
obj = gmdistribution.fit(imgvect,2,'Options',gmopt);
% Only I want to go unimodal if bimodal the exited on warning
if strcmp(lastwarn,['Failed to converge in ',num2str(numIter),' iterations for gmdistribution with 2 components'])
    obj = gmdistribution.fit(imgvect,1,'Options',gmopt);     
end    
lastwarn('','TestEnv:InvalidInput'); % Reset warning to something else random  
% Get the minimum mean and its standard deviation (background) and normalize that to [0,1]
mu_ind = find(obj.mu == min(obj.mu),1,'first');
background_mu = obj.mu(mu_ind);
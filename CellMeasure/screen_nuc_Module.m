function [CellMeasurements, ModuleData] = screen_nuc_Module(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = screen_nuc_Module(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCREEN1MODULE is closely related to tracking-based nuclear measurement modules, but is tweaked slightly to account
% for screening (e.g. applies universal correction to image background)
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;
measure_cc = label2cc(labels.Nucleus,0);

% Mode-balance 1st auxililiary image - bimodal distribution assumed (nuclear expression, cytoplasmic expression, and b.g.)
corr_img = AuxImages{1};
% Background correct
if isfield(ModuleData,'X')
    warning off MATLAB:nearlySingularMatrix
    pStar = (X'*X)\(X')*corr_img(:);
    % Apply correction
    corr_img = reshape((double(corr_img(:) - X*pStar)),size(corr_img));
end

corr_img(imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius))) = []; % Drop foreground objects for correction calculation
[~, dist1] = modebalance(corr_img,2,ModuleData.BitDepth,'measure'); 
AuxImages{1} = (AuxImages{1} - dist1(1))/dist1(2); % Background subtract/divide

% Intensity-based measurement initialization
CellMeasurements.ImageBackground1 = dist1';
CellMeasurements.MeanNuc1 =  nan(parameters.TotalCells,parameters.TotalImages);
CellMeasurements.IntegratedNuc1 =  nan(parameters.TotalCells,parameters.TotalImages);
CellMeasurements.MedianNuc1 = nan(parameters.TotalCells,parameters.TotalImages);


% Cycle through each cell and assign measurements
for n = 1:measure_cc.NumObjects
    CellMeasurements.MeanNuc1(n,iteration) = nanmean(AuxImages{1}(measure_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedNuc1(n,iteration) = nansum(AuxImages{1}(measure_cc.PixelIdxList{n}));
    CellMeasurements.MedianNuc1(n,iteration) = nanmedian(AuxImages{1}(measure_cc.PixelIdxList{n}));
end


% Measure nuclei in 2nd auxiliary, if it is specified
if ~isempty(AuxImages{2})
    % Mode-balance 1st auxililiary image - bimodal distribution assumed (nuclear expression, cytoplasmic expression, and b.g.)
    corr_img = AuxImages{2};
    % Background correct
    if isfield(ModuleData,'X')
        warning off MATLAB:nearlySingularMatrix
        pStar = (X'*X)\(X')*corr_img(:);
        % Apply correction
        corr_img = reshape((double(corr_img(:) - X*pStar)),size(corr_img));
    end
    corr_img(imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius))) = []; % Drop foreground objects for correction calculation
    [~, dist1] = modebalance(corr_img,2,ModuleData.BitDepth,'measure'); 
    AuxImages{2} = (AuxImages{2} - dist1(1))/dist1(2); % Background subtract/divide

    % Intensity-based measurement initialization
    CellMeasurements.ImageBackground2 = dist1';
    CellMeasurements.MeanNuc2 =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedNuc2 =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianNuc2 = nan(parameters.TotalCells,parameters.TotalImages);

    % Cycle through each cell and assign measurements
    for n = 1:measure_cc.NumObjects
        CellMeasurements.MeanNuc2(n,iteration) = nanmean(AuxImages{2}(measure_cc.PixelIdxList{n}));
        CellMeasurements.IntegratedNuc2(n,iteration) = nansum(AuxImages{2}(measure_cc.PixelIdxList{n}));
        CellMeasurements.MedianNuc2(n,iteration) = nanmedian(AuxImages{2}(measure_cc.PixelIdxList{n}));
    end

end
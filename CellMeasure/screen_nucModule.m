function [CellMeasurements, ModuleData] = screen_nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = screen_nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCREEN_NUCMODULE is a fixed-cell measurement module. It will perform flatfield/background correction for each
% specified image (Measure img 1 corresponds to Flatfield{1}, measure image 2 with Flatfield{2}, etc.), then will
% calculate nuclear expression levels for each individual cell.
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

measure_cc = label2cc(labels.Nucleus,0);


for i = 1:length(AuxImages)
    if ~isempty(AuxImages{i})
        % STEP 1: image normalization
        % Normalization 1: flatfield correction -> estimate range for scaling
        corr_img = flatfieldcorrect(double(AuxImages{i}),double(parameters.Flatfield{i}));
        % Normalization 2: mode-balance - bimodal distribution assumed after dropping nuclei (leaves cytoplasmic + b.g.)
        corr_img = corr_img - min(corr_img(:));
        tmp = corr_img;
        tmp(imdilate(labels.Nucleus>0,diskstrel(round(parameters.MinNucleusRadius*2)))) = []; % Drop foreground objects for correction calculation
        [~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure');
        corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
        AuxImages{i} = corr_img;
        
        % STEP 2: initialize data
        CellMeasurements.(['MeanNuc',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['IntegratedNuc',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['MedianNuc',num2str(i)]) = nan(parameters.TotalCells,parameters.TotalImages,1);

        % Step 3: measure data
        for n = 1:measure_cc.NumObjects
            CellMeasurements.(['MeanNuc',num2str(i)])(n) = nanmean(AuxImages{i}(measure_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedNuc',num2str(i)])(n) = nansum(AuxImages{i}(measure_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianNuc',num2str(i)])(n) = nanmedian(AuxImages{i}(measure_cc.PixelIdxList{n}));
        end
    end
end
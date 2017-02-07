function [CellMeasurements, ModuleData] = screen_cytoModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = screen_cytoModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCREEN_NUCMODULE is a fixed-cell measurement module. It will perform flatfield/background correction for each
% specified image (Measure img 1 corresponds to Flatfield{1}, measure image 2 with Flatfield{2}, etc.), then will
% calculate cytoplasmic and whole-cell expression levels for each individual cell.
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% If cells were not segmented; do annulus to estimate cytoplasmic area.
if isequal(labels.Cell,labels.Nucleus)
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*2));
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end

% Make conncomp structures
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:cyto_cc.NumObjects
    cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
end


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
        CellMeasurements.(['MeanCell',num2str(i)]) =  nan(parameters.TotalCells,1);
        CellMeasurements.(['IntegratedCell',num2str(i)]) =  nan(parameters.TotalCells,1);
        CellMeasurements.(['MedianCell',num2str(i)]) = nan(parameters.TotalCells,1);
        CellMeasurements.(['MeanCyto',num2str(i)]) =  nan(parameters.TotalCells,1);
        CellMeasurements.(['IntegratedCyto',num2str(i)]) =  nan(parameters.TotalCells,1);
        CellMeasurements.(['MedianCyto',num2str(i)]) = nan(parameters.TotalCells,1);
        
        
        % Step 3: measure data
        for n = 1:cell_cc.NumObjects
            CellMeasurements.(['MeanCell',num2str(i)])(n) = nanmean(AuxImages{i}(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedCell',num2str(i)])(n) = nansum(AuxImages{i}(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianCell',num2str(i)])(n) = nanmedian(AuxImages{i}(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['MeanCyto',num2str(i)])(n) = nanmean(AuxImages{i}(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedCyto',num2str(i)])(n) = nansum(AuxImages{i}(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianCyto',num2str(i)])(n) = nanmedian(AuxImages{i}(cyto_cc.PixelIdxList{n}));
        end
    end
end
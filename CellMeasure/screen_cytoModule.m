function [CellMeasurements, ModuleData] = screen_cytoModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = screen_cytoModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCREEN_CYTOMODULE is a fixed-cell measurement module. It will perform flatfield/background correction for each
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
if strcmpi(parameters.ImageType,'none')
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*2));
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end

% Make conncomp structures
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:cyto_cc.NumObjects
    if ~isempty(cyto_cc.PixelIdxList{i})
        cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
    end
end


for i = 1:length(AuxImages)
    if ~isempty(AuxImages{i})
        % STEP 1: image normalization
        % Normalization 1: flatfield correction -> estimate range for scaling
        corr_img = flatfieldcorrect(double(AuxImages{i}),double(parameters.Flatfield{i}));
        % Normalization 2: mode-balance - unimodal distribution assumed after dropping cells (b.g. only)
        corr_img = corr_img - prctile(corr_img(:),1); % Need to have some kind of confluent "switch" here?
%         tmp = corr_img;
%         tmp(labels.Cell>0) = []; % Drop foreground objects for correction calculation
%         if ~isempty(tmp)
%             [~, dist1] = modebalance(tmp,0,ModuleData.BitDepth,'measure'); % Get mode of existing background
%         else
%             dist1 = prctile(corr_img(:),2); % If there is no background (100% confluence), just subtract 2nd percentile
%         end
%         corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
        AuxImages{i} = corr_img;
         % STEP 2: initialize data
        CellMeasurements.(['MeanCell',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['IntegratedCell',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['MedianCell',num2str(i)]) = nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['MeanCyto',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['IntegratedCyto',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
        CellMeasurements.(['MedianCyto',num2str(i)]) = nan(parameters.TotalCells,parameters.TotalImages,1);

        % STEP 3: measure data
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
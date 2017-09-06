function [CellMeasurements, ModuleData] = ktrModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% KTRMODULE  measure cytoplasmic:nuclear ratio in auxiliary (fluorescent) image
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% If cells were not segmented, do annulus to estimate cytoplasmic area.
if strcmpi(parameters.ImageType,'none')
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*2));
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end

% Make conncomp structures for nucleus, cell, and cytoplasm
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:cyto_cc.NumObjects
    if ~isempty(cyto_cc.PixelIdxList{i})
        cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
    end
end

% Flatfield-correct auxiliary image
for i = 1:length(AuxImages)
    if ~isempty(AuxImages{i})
        
        % STEP 1: IMAGE NORMALIZATION
        % Normalization 1: flatfield correction -> estimate range for scaling
        corr_img = flatfieldcorrect(double(AuxImages{i})-double(parameters.Flatfield{end}),double(parameters.Flatfield{i}));
        % Normalization 2: background-subtract -> assume image is pretty confluent...
        corr_img = corr_img - prctile(corr_img(:),2);        
        AuxImages{i} = corr_img;
        
        
        % STEP 2: INITIALIZE DATA
        if ~isfield(CellMeasurements,['KTR_ratio',num2str(i)])
            CellMeasurements.(['KTR_ratio',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
            CellMeasurements.(['KTR_nuc',num2str(i)]) =  nan(parameters.TotalCells,parameters.TotalImages,1);
            CellMeasurements.(['KTR_cyto',num2str(i)]) = nan(parameters.TotalCells,parameters.TotalImages,1);
        end

        % STEP 3: MAKE MEASUREMENT
        for n = 1:cyto_cc.NumObjects        
            if ~isempty(nuc_cc.PixelIdxList{n})
                % (Median of nuclear and cytoplasmic compartments)
                median_n = median(AuxImages{i}(nuc_cc.PixelIdxList{n}));
                median_c = median(AuxImages{i}(cyto_cc.PixelIdxList{n}));
                % Assign measurements
                CellMeasurements.(['KTR_nuc',num2str(i)])(n,ModuleData.iter) = median_n;
                CellMeasurements.(['KTR_cyto',num2str(i)])(n,ModuleData.iter) = median_c;
                CellMeasurements.(['KTR_ratio',num2str(i)])(n,ModuleData.iter) = median_c/median_n;
            end 
        end 
    end
end


%% OLD normalization (background modeling and subtraction)
%         tmp = corr_img;
%         tmp(imdilate(labels.Nucleus>0,diskstrel(round(parameters.MinNucleusRadius*2)))) = []; % Drop foreground objects for correction calculation
%         [~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure');
%         corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)

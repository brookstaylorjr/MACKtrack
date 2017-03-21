function [CellMeasurements, ModuleData] = intensityModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INTENSITYMODULE measures cellular intensity (nuclear, cell, and cytoplasmic) in specified AuxImages channels
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Get indicies and bwconncomp structures for measurement
iteration  = ModuleData.iter;

% Make conncomp structures for nucleus, cell, and cytoplasm
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:cyto_cc.NumObjects
    if ~isempty(cyto_cc.PixelIdxList{i})
        cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
    end
end



% Main cycle: correct image, initialize data (if not present), make measurments
for img = 1:length(AuxImages)
    if ~isempty(AuxImages{img})
        % 1) Background correct image (try to do flatfield, if available)
        if (length(parameters.Flatfield)>=img) && isequal(size(AuxImages{img}),size(parameters.Flatfield{img}))
            img0 = flatfieldcorrect(AuxImages{img},double(parameters.Flatfield{img}));
            img0 = img0-prctile(img0(:),2); % Background subtract
        else
            if ~isfield(ModuleData,'distr')
                [img0, ModuleData.distr] = modebalance(AuxImages{img},2,ModuleData.BitDepth,'measure'); 
            else
                img0 = modebalance(AuxImages{img},2,ModuleData.BitDepth,'correct',ModuleData.distr);
            end
        end


        % - - - - NUCLEAR measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,['MeanIntensity_nuc',num2str(img)])
            CellMeasurements.(['MeanIntensity_nuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['MedianIntensity_nuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['IntegratedIntensity_nuc',num2str(img)]) = nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:nuc_cc.NumObjects
            CellMeasurements.(['MeanIntensity_nuc',num2str(img)])(n,iteration) = nanmean(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianIntensity_nuc',num2str(img)])(n,iteration) = nanmedian(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedIntensity_nuc',num2str(img)])(n,iteration) = nansum(img0(cyto_cc.PixelIdxList{n}));
        end



        % - - - - CYTOPLASMIC/WHOLE-CELL measurements - - - -
        if ~strcmpi(parameters.ImageType,'none')
            % A) Initialize fields
            if ~isfield(CellMeasurements,(['MeanIntensity_cyto',num2str(img)]))
                CellMeasurements.(['MeanIntensity_cyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
                CellMeasurements.(['MedianIntensity_cyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
                CellMeasurements.(['IntegratedIntensity_cyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);

                CellMeasurements.(['MeanIntensity_cell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
                CellMeasurements.(['MedianIntensity_cell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
                CellMeasurements.(['IntegratedIntensity_cell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            end

            % B) Assign measurements
            for n = 1:cyto_cc.NumObjects
                CellMeasurements.(['MeanIntensity_cyto',num2str(img)])(n,iteration) = nanmean(img0(cyto_cc.PixelIdxList{n}));
                CellMeasurements.(['MedianIntensity_cyto',num2str(img)])(n,iteration) = nansum(img0(cyto_cc.PixelIdxList{n}));
                CellMeasurements.(['IntegratedIntensity_cyto',num2str(img)])(n,iteration) = nanmedian(img0(cyto_cc.PixelIdxList{n}));

                CellMeasurements.(['MeanIntensity_cell',num2str(img)])(n,iteration) = nanmean(img0(cell_cc.PixelIdxList{n}));
                CellMeasurements.(['MedianIntensity_cell',num2str(img)])(n,iteration) = nansum(img0(cell_cc.PixelIdxList{n}));
                CellMeasurements.(['IntegratedIntensity_cell',num2str(img)])(n,iteration) = nanmedian(img0(cell_cc.PixelIdxList{n}));
            end
        end
    end
end




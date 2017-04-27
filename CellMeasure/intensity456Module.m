function [CellMeasurements, ModuleData] = intensity456Module(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INTENSITY456MODULE measures cellular intensity (nuclear, cell, and cytoplasmic) in specified AuxImages channels
% (Allows extra-channel measurements, e.g. MeanCell4, MeanCell5, and MeanCell6) -> normalization will use Flatfield4
% Flatfield5, and Flatfield6, accordingly
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
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*1.5)); % (Expand by 1.5 radii)
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end


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
        if (length(parameters.Flatfield)>=(img+3)) && isequal(size(AuxImages{img}),size(parameters.Flatfield{img+3}))
            img0 = flatfieldcorrect(AuxImages{img},double(parameters.Flatfield{img+3}));
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
        if ~isfield(CellMeasurements,['MeanNuc',num2str(img+3)])
            CellMeasurements.(['MeanNuc',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['MedianNuc',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['IntegratedNuc',num2str(img+3)]) = nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:nuc_cc.NumObjects
            CellMeasurements.(['MeanNuc',num2str(img+3)])(n,iteration) = nanmean(img0(nuc_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianNuc',num2str(img+3)])(n,iteration) = nanmedian(img0(nuc_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedNuc',num2str(img+3)])(n,iteration) = nansum(img0(nuc_cc.PixelIdxList{n}));
        end



        % - - - - CYTOPLASMIC/WHOLE-CELL measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,(['MeanCyto',num2str(img+3)]))
            CellMeasurements.(['MeanCyto',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['MedianCyto',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['IntegratedCyto',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);

            CellMeasurements.(['MeanCell',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['MedianCell',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['IntegratedCell',num2str(img+3)]) =  nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:cyto_cc.NumObjects
            CellMeasurements.(['MeanCyto',num2str(img+3)])(n,iteration) = nanmean(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianCyto',num2str(img+3)])(n,iteration) = nanmedian(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedCyto',num2str(img+3)])(n,iteration) = nansum(img0(cyto_cc.PixelIdxList{n}));

            CellMeasurements.(['MeanCell',num2str(img+3)])(n,iteration) = nanmean(img0(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['MedianCell',num2str(img+3)])(n,iteration) = nanmedian(img0(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['IntegratedCell',num2str(img+3)])(n,iteration) = nansum(img0(cell_cc.PixelIdxList{n}));
        end 
    end
end




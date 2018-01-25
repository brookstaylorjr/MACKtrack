function [CellMeasurements, ModuleData] = cellcycleModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CELLCYCLE measures cellular intensity (nuclear, cell, and cytoplasmic) in specified AuxImages channels - 
% written similarly to intensityModule, but with:
%      - more aggressive background normalization (subtracts mean of unimodal distr. fitted to BG pixels)
%      - stricter cell boundaries (thinner annulus)
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
    tmp_mask = imdilate(labels.Nucleus>0,ones(5)); % (Expand by 5 pixels in all directions)
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

% Define background region
drop_region = imdilate(labels.Nucleus>0,diskstrel(2*parameters.MinNucleusRadius));

% Main cycle: correct image, initialize data (if not present), make measurments
for img = 1:length(AuxImages)
    if ~isempty(AuxImages{img})
        % 1) Background correct image (try to do flatfield, if available)
        if (length(parameters.Flatfield)>img)
            if isequal(size(AuxImages{img}),size(parameters.Flatfield{img}))
                AuxImages{img} = double(AuxImages{img}) - double(parameters.Flatfield{end});
                img0 = flatfieldcorrect(AuxImages{img},double(parameters.Flatfield{img}));
                % Fit unimodal normal distribution to image background. Subtract that mean from original image
                [~, d] = modebalance(img0(~drop_region),1,ModuleData.BitDepth,'measure');
                img0 = img0-d(1);
            else
                error(['Size mismatch between provided flatfield (#', num2str(img), ' and AuxImage'])
            end
        else
            if ~isfield(ModuleData,'distr')
                [img0, ModuleData.distr] = modebalance(AuxImages{img},2,ModuleData.BitDepth,'measure'); 
            else
                img0 = modebalance(AuxImages{img},2,ModuleData.BitDepth,'correct',ModuleData.distr);
            end
        end


        % - - - - NUCLEAR measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,['CCMeanNuc',num2str(img)])
            CellMeasurements.(['CCMeanNuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCMedianNuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCIntegratedNuc',num2str(img)]) = nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:nuc_cc.NumObjects
            CellMeasurements.(['CCMeanNuc',num2str(img)])(n,iteration) = nanmean(img0(nuc_cc.PixelIdxList{n}));
            CellMeasurements.(['CCMedianNuc',num2str(img)])(n,iteration) = nanmedian(img0(nuc_cc.PixelIdxList{n}));
            CellMeasurements.(['CCIntegratedNuc',num2str(img)])(n,iteration) = nansum(img0(nuc_cc.PixelIdxList{n}));
        end



        % - - - - CYTOPLASMIC/WHOLE-CELL measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,(['CCMeanCyto',num2str(img)]))
            CellMeasurements.(['CCMeanCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCMedianCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCIntegratedCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);

            CellMeasurements.(['CCMeanCell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCMedianCell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['CCIntegratedCell',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:cyto_cc.NumObjects
            CellMeasurements.(['CCMeanCyto',num2str(img)])(n,iteration) = nanmean(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['CCMedianCyto',num2str(img)])(n,iteration) = nanmedian(img0(cyto_cc.PixelIdxList{n}));
            CellMeasurements.(['CCIntegratedCyto',num2str(img)])(n,iteration) = nansum(img0(cyto_cc.PixelIdxList{n}));

            CellMeasurements.(['CCMeanCell',num2str(img)])(n,iteration) = nanmean(img0(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['CCMedianCell',num2str(img)])(n,iteration) = nanmedian(img0(cell_cc.PixelIdxList{n}));
            CellMeasurements.(['CCIntegratedCell',num2str(img)])(n,iteration) = nansum(img0(cell_cc.PixelIdxList{n}));
        end
        
    end
end





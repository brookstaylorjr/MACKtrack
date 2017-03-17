function [CellMeasurements, ModuleData] = intensityModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INTENSITYMODULE measures cellular intensity in AuxImages{1} channel
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
nuc_cc = label2cc(labels.Nucleus,0);

if isfield(labels,'Cell')
    cyto_label = labels.Cell;
    cyto_label(labels.Nucleus>0) = 0; % cytoplasm only
    cyto_cc = label2cc(cyto_label,0);
    cell_cc = label2cc(labels.Cell,0);
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
        % 
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
    if isfield(labels,'Cell')
        % A) Initialize fields
        if ~isfield(CellMeasurements,(['MeanIntensity_cyto',num2str(img)]))
            CellMeasurements.(['MedianIntensity_cyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['IntegratedIntensity_cyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);

            CellMeasurements.MeanFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.IntegratedFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.MedianFRET_cell = nan(parameters.TotalCells,parameters.TotalImages);

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





% Cycle through each cell and assign measurements
cells = unique(labels.Cell(labels.Cell>0));
for n = 1:length(cells)
    CellMeasurements.MeanIntensity(cells(n),iteration) = mean(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.IntegratedIntensity(cells(n),iteration) = sum(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.MedianIntensity(cells(n),iteration) = median(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.IntensityPercentiles(cells(n),iteration,:) = prctile(AuxImages{1}(labels.Cell==cells(n)),1:100);
end


% Measure cells in 2nd auxiliary image, if it is specified
if ~isempty(AuxImages{2})
    % Mode-balance
    if ~isfield(ModuleData,'distr2')
        [~, ModuleData.distr2] = modebalance(AuxImages{2},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{2} = modebalance(AuxImages{2},2,ModuleData.BitDepth,'correct',ModuleData.distr2);
    end


    % On first call, initialize all new CellMeasurements fields 
    if ~isfield(CellMeasurements,'MeanIntensityNuc2')
        % Intensity-based measurement initialization
        CellMeasurements.MeanIntensity2 =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.IntegratedIntensity2 =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.MedianIntensity2 = nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.IntensityPercentiles2 = nan(parameters.TotalCells,parameters.TotalImages,100);

    end

    % Cycle through each cell and assign measurements
    cells = unique(labels.Cell(labels.Cell>0));
    for n = 1:length(cells)
        CellMeasurements.MeanIntensity2(cells(n),iteration) = mean(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.IntegratedIntensity2(cells(n),iteration) = sum(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.MedianIntensity2(cells(n),iteration) = median(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.IntensityPercentiles2(cells(n),iteration,:) = prctile(AuxImages{2}(labels.Cell==cells(n)),1:100);

    end

end

ModuleDataOut = ModuleData;
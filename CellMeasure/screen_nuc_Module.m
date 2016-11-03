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
fun = @(block_struct) prctile(block_struct.data(:),2);

% Background correct - find best fit scale from 0x to 10x, 0.5x increment
stdev1 = inf;
for i = 0:20
    bg_subtract = AuxImages{1}-ModuleData.Flatfield{1}*0.5*i;
    lo_vals = blockproc(bg_subtract,[200 200],fun);
    if std(lo_vals(:))<stdev1
        stdev1 = std(lo_vals(:));
        mult = 0.5*i;
    end
end
% Mode-balance image - bimodal distribution assumed after dropping foreground objects (cytoplasmic expression, and b.g.)
corr_img = AuxImages{1} - ModuleData.Flatfield{1}*mult;
corr_img = corr_img - min(corr_img(:));
tmp = corr_img;
tmp(imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*4))) = []; % Drop foreground objects for correction calculation
[~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure'); 
corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
AuxImages{1} = corr_img;

% Intensity-based measurement initialization
CellMeasurements.Mult1 = mult;
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
    stdev1 = inf;
    for i = 0:20
        bg_subtract = AuxImages{2}-ModuleData.Flatfield{1}*0.25*i;
        lo_vals = blockproc(bg_subtract,[200 200],fun);
        if std(lo_vals(:))<stdev1
            stdev1 = std(lo_vals(:));
            mult = 0.25*i;
        end
    end
    corr_img = AuxImages{2} - ModuleData.Flatfield{1}*mult;
    corr_img = corr_img - min(corr_img(:));
    tmp = corr_img;
    tmp(imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*4))) = []; % Drop foreground objects for correction calculation
    [~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure'); 
    corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
    AuxImages{2} = corr_img;

    % Intensity-based measurement initialization
    CellMeasurements.Mult2 = mult;
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
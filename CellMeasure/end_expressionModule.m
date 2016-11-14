function [CellMeasurements, ModuleData] = end_expressionModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = end_expressionModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% END_EXPRESSIONMODULE performs endopoint expression measurement (e.g. post staining of the same cells being imaged)
% This involves (1) making sure cell wasn't washed away, (2) performing segmentation on AuxImages{1}, and (3) making
% basic whole-cell expression measurements in AuxImages{1} and {2}. ONLY triggered @ last frame.
%
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ModuleData.iter == parameters.TotalImages
    %% Load the corresponding (final) nuclear image - calculate the "jump" to the endpoint images
    home_folder = mfilename('fullpath');
    slash_idx = strfind(home_folder,filesep);
    load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
    i = ModuleData.i; j = ModuleData.j;
    nucName1 = eval(parameters.NucleusExpr);
    images.nuc = checkread([locations.scope,parameters.ImagePath,nucName1],ModuleData.BitDepth,1,parameters.debug);
    jump1 = calculatejump(images.nuc,AuxImages{1});
    
    % Translate corresponding nuclear label matrix
	nuc_label = imtranslate(labels.Nucleus,[jump1(2),jump1(1)]);
    
    % Flatfield-correct aux image #2 and make cell mask
    cell_img = flatfieldcorrect(double(AuxImages{2}),double(parameters.Flatfield{1}));
    cell_img = cell_img-min(cell_img(:))+2;
    mask0 = cell_img>tsaithresh(cell_img,false(size(cell_img)));
    % Morphological cleanup on mask
    mask1 = bwareaopen(mask0,2);
    mask1 = imclose(mask1,diskstrel(3));
    mask1 = bwareaopen(mask1,parameters.NoiseSize);
    mask1 = ~bwareaopen(~mask1,parameters.MinHoleSize);

    % Zero out any nuclei that got washed away (i.e nucleus no longer overlaps with a cell)
    nuc_ids = unique(nuc_label(nuc_label>0));
    overlap = zeros(max(nuc_ids),1);
    tmp_cc = label2cc(nuc_label,0);
    for n = 1:length(nuc_ids)
        if ~isempty(tmp_cc.PixelIdxList{nuc_ids(n)})
            overlap(nuc_ids(n)) = sum(mask1(tmp_cc.PixelIdxList{nuc_ids(n)}))./numel(tmp_cc.PixelIdxList{nuc_ids(n)});
        end
    end    
    nuc_ids = find(overlap>0.66);
    nuc_label(~ismember(nuc_label,nuc_ids)) = 0; 
    cell_label = propagatesegment(nuc_label,mask1,AuxImages{1},2);
    
    % Make a diagnostic output
    save_dir = [locations.data,filesep,parameters.SaveDirectory,filesep,'EndpointSegmentation',filesep];
    if ~exist(save_dir,'dir');  mkdir(save_dir); end
   saveFig(cell_img,cell_label,nuc_label,[],ModuleData.BitDepth,...
        [save_dir,'Endpoint_pos',numseq(ModuleData.i,2),'.jpg'],0.2,[1024 1024], [-3 60]);
    
    % Background subtraction on cell image
    [~, dist1] = modebalance(cell_img,2,ModuleData.BitDepth,'measure'); 
    cell_img = (cell_img - dist1(1)); % Background subtract (DON'T divide) 
    
    % Initialize measurements for cell_img (AuxImages{2})
    CellMeasurements.EndMean1 =  nan(parameters.TotalCells,1);
    CellMeasurements.EndIntegrated1 =  nan(parameters.TotalCells,1);
    CellMeasurements.EndMedian1 = nan(parameters.TotalCells,1);

    % Cycle objects and measure
    measure_cc = label2cc(cell_label,0);
    for n = 1:measure_cc.NumObjects
        CellMeasurements.EndMean1(n) = nanmean(cell_img(measure_cc.PixelIdxList{n}));
        CellMeasurements.EndIntegrated1(n) = nansum(cell_img(measure_cc.PixelIdxList{n}));
        CellMeasurements.EndMedian1(n) = nanmedian(cell_img(measure_cc.PixelIdxList{n}));
    end
    
    if ~isempty(AuxImages{3})
        % Flatfield/background subtraction correction
        cell_img2 = flatfieldcorrect(double(AuxImages{3}),double(parameters.Flatfield{1}));
        cell_img2 = cell_img2-min(cell_img2(:))+2;
        [~, dist1] = modebalance(cell_img2,2,ModuleData.BitDepth,'measure'); 
        cell_img2 = (cell_img2 - dist1(1)); % Background subtract (DON'T divide) 
        
         % Initialize measurements for AuxImages{2}
        CellMeasurements.EndMean2 =  nan(parameters.TotalCells,1);
        CellMeasurements.EndIntegrated2 =  nan(parameters.TotalCells,1);
        CellMeasurements.EndMedian2 = nan(parameters.TotalCells,1);

        % Cycle objects and measure
        measure_cc = label2cc(cell_label,0);
        for n = 1:measure_cc.NumObjects
            CellMeasurements.EndMean2(n) = nanmean(cell_img2(measure_cc.PixelIdxList{n}));
            CellMeasurements.EndIntegrated2(n) = nansum(cell_img2(measure_cc.PixelIdxList{n}));
            CellMeasurements.EndMedian2(n) = nanmedian(cell_img2(measure_cc.PixelIdxList{n}));
        end
    end
    



end
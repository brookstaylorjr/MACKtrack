function [CellMeasurements, ModuleData] = textureModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% textureMODULE measures texture and contrast of specified AuxImages channels
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

% create mask so we can make background nan;
nuc_mask = labels.Nucleus>0;
cell_mask = labels.Cell>0;
cyto_mask = cell_mask - nuc_mask;

% Get indicies and bwconncomp structures for measurement
iteration  = ModuleData.iter;

% Make conncomp structures for nucleus, cell, and cytoplasm
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:min([length(cyto_cc.PixelIdxList),length(nuc_cc.PixelIdxList)])
    if ~isempty(cyto_cc.PixelIdxList{i})
        cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
    end
end

% Main cycle: correct image, initialize data (if not present), make measurments
for img = 1:length(AuxImages)
    if ~isempty(AuxImages{img})
        % 1) Background correct image (try to do flatfield, if available)
        if (length(parameters.Flatfield)>img)
            if isequal(size(AuxImages{img}),size(parameters.Flatfield{img}))
                AuxImages{img} = double(AuxImages{img}) - double(parameters.Flatfield{end});
                img0 = flatfieldcorrect(AuxImages{img},double(parameters.Flatfield{img}));
                img0 = img0-prctile(img0(:),2); % Background subtract
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
        
        nucimg = img0;
        nucimg(~nuc_mask) = nan;
        cytoimg = img0;
        cytoimg(~cyto_mask) = nan;

        % - - - - NUCLEAR measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,['TextureNuc',num2str(img)])
            CellMeasurements.(['TextureNuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['ContrastNuc',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['StdevNuc',num2str(img)]) = nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        offset1 = [0 1; -1 1; -1 0;-1 -1]; % offset represents spatial feature measurements in all directions;
        for n = 1:nuc_cc.NumObjects
            [imrow,imcol] = sub2ind(size(img0),[nuc_cc.PixelIdxList{n}]);
            cropimg0nuc = nucimg(min(imrow):max(imrow),min(imcol):max(imcol));
            [graycomat,~] = graycomatrix(cropimg0nuc,'Offset',offset1,'GrayLimits',[min(cropimg0nuc(:)) max(cropimg0nuc(:))],'NumLevels',64,'Symmetric',true);
            graystats = graycoprops(graycomat,{'correlation','contrast'});
            CellMeasurements.(['TextureNuc',num2str(img)])(n,iteration) = nanmean(graystats.Correlation);
            CellMeasurements.(['ContrastNuc',num2str(img)])(n,iteration) = nanmean(graystats.Contrast);
            CellMeasurements.(['StdevNuc',num2str(img)])(n,iteration) = nanstd(img0(nuc_cc.PixelIdxList{n}));
        end



        % - - - - CYTOPLASMIC/WHOLE-CELL measurements - - - -
        % A) Initialize fields
        if ~isfield(CellMeasurements,(['TextureCyto',num2str(img)]))
            CellMeasurements.(['TextureCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['ContrastCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
            CellMeasurements.(['StdevCyto',num2str(img)]) =  nan(parameters.TotalCells,parameters.TotalImages);
        end

        % B) Assign measurements
        for n = 1:cyto_cc.NumObjects
            [imrow,imcol] = sub2ind(size(img0),[cyto_cc.PixelIdxList{n}]);
            cropimg0cyto = cytoimg(min(imrow):max(imrow),min(imcol):max(imcol));
            [graycomat,~] = graycomatrix(cropimg0cyto,'Offset',offset1,'GrayLimits',[min(cropimg0cyto(:)) max(cropimg0cyto(:))],'NumLevels',64,'Symmetric',true);
            graystats = graycoprops(graycomat,{'correlation','contrast'});
            CellMeasurements.(['TextureCyto',num2str(img)])(n,iteration) = nanmean(graystats.Correlation);
            CellMeasurements.(['ContrastCyto',num2str(img)])(n,iteration) = nanmean(graystats.Contrast);
            CellMeasurements.(['StdevCyto',num2str(img)])(n,iteration) = nanstd(img0(nuc_cc.PixelIdxList{n}));
        end
        
    end
end





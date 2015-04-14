function [CellMeasurements, ModuleDataOut] = nfkbnucModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMODULE  measure single-cell nuclear fraction in auxiliary image
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Mode-balance 1st auxililiary image - bimodal distribution assumed
    if ~isfield(ModuleData,'distr')
        [AuxImages{1}, ModuleData.distr] = modebalance(AuxImages{1},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{1} = modebalance(AuxImages{1},2,ModuleData.BitDepth,'correct',ModuleData.distr);
    end
nfkb = AuxImages{1};


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'NFkBFraction')
    % Intensity-based measurement initialization
    CellMeasurements.NFkBNuclear = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.NFkBCytoplasm = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.NFkBFraction =  zeros(parameters.TotalCells,parameters.TotalImages);
end

% Count cells for each frame, initialize bins
cells = unique(labels.Nucleus);
cells(cells==0) = [];

% Cytoplasm wasn't segmented, so we need to "form" it
labels.Cell = imdilate(labels.Nucleus,diskstrel(10));
labels.Cell(labels.Nucleus>0) = labels.Nucleus(labels.Nucleus>0);



% Cycle through each image and assign measurements
for i = 1:length(cells)
    bins = 0:255;
    nucleus = imerode(labels.Nucleus==cells(i),diskstrel(5));
    cytoplasm = (labels.Cell==cells(i)) &~nucleus;
    
    % Use median of nuclei compartment (mode would probably also work just fine)
    median_n = median(nfkb(nucleus));

    % Get secondary (higher) mode of cytoplasmic compartment
    n = hist(nfkb(cytoplasm),bins);
    thresh1 = otsuthresh(nfkb,~cytoplasm,'none');
    if ~isempty(thresh1)
        n(bins<=thresh1) = [];
        bins(bins<=thresh1) = [];
    end
    if ~isempty(n)
        mode2_c = bins(n==max(n));
    else
        mode2_c = median(nfkb(cytoplasm));
    end

    % Assign measurements
    CellMeasurements.NFkBNuclear(cells(i),ModuleData.iter) = median_n(1);
    CellMeasurements.NFkBCytoplasm(cells(i),ModuleData.iter) = mode2_c(1);
    CellMeasurements.NFkBFraction(cells(i),ModuleData.iter) = median_n(1)/mode2_c(1);
end

% On last iteration, normalize nuclear measurement using 10th percentile values 
if ModuleData.iter == parameters.TotalImages
    % Get each cell's 10th percentile
    bottoms = prctile(CellMeasurements.NFkBNuclear,10,2);
    for i = 1:size(CellMeasurements.NFkBNuclear,1)
        CellMeasurements.NFkBNuclear(i,:) = CellMeasurements.NFkBNuclear(i,:) - bottoms(i);
    end
end

ModuleDataOut = ModuleData;

function [CellMeasurements, ModuleDataOut] = nfkbdimModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBDIMMODULE  measure nuclear fraction in auxiliary (fluorescent) image.
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Get NFkB image and background-correct
nfkb = AuxImages{1};
if ~isfield(ModuleData,'X')
    ModuleData.X = backgroundcalculate(size(nfkb));
end

warning off MATLAB:nearlySingularMatrix
pStar = (ModuleData.X'*ModuleData.X)\(ModuleData.X')*double(nfkb(:));
warning on MATLAB:nearlySingularMatrix

% Apply background correction
nfkb = reshape((double(nfkb(:) - ModuleData.X*pStar)),size(nfkb));
nfkb = nfkb-min(nfkb(:)); % Set minimum to zero

if ~isfield(ModuleData,'distr')
    [nfkb, ModuleData.distr] = modebalance(nfkb,1,ModuleData.BitDepth,'measure');
else
    nfkb = modebalance(nfkb,1,ModuleData.BitDepth,'correct',ModuleData.distr);
end


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'NFkBdimNuclear')
    % Intensity-based measurement initialization
    CellMeasurements.NFkBdimNuclear = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.NFkBdimCytoplasm = nan(parameters.TotalCells,parameters.TotalImages);
end

% Count cells for each frame, initialize bins
cells = unique(labels.Nucleus);
cells(cells==0) = [];

% Cycle through each image and assign measurements
rng = [0 max(nfkb(:))];
bins = rng(1):round(diff(rng)/255):rng(2);
for i = 1:length(cells)
    nucleus = imerode(labels.Nucleus==cells(i),diskstrel(floor(parameters.MinNucleusRadius/5)));
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
    CellMeasurements.NFkBdimNuclear(cells(i),ModuleData.iter) = median_n(1);
    CellMeasurements.NFkBdimCytoplasm(cells(i),ModuleData.iter) = mode2_c(1);
end

ModuleDataOut = ModuleData;

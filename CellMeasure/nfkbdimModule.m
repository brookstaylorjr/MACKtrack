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
    [~, ModuleData.distr] = modebalance(nfkb,1,ModuleData.BitDepth,'measure');
else
    nfkb = modebalance(nfkb,1,ModuleData.BitDepth,'correct',ModuleData.distr);
end

% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'NFkBdimNuclear')
    % Intensity-based measurement initialization
    CellMeasurements.NFkBdimNuclear = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.NFkBdimNuclear_erode = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.NFkBdimCytoplasm = nan(parameters.TotalCells,parameters.TotalImages);
end

% Count cells for each frame, initialize bins
cells = unique(labels.Nucleus);
cells(cells==0) = [];

[~, distr2] = modebalance(nfkb(labels.Cell==0),1,ModuleData.BitDepth,'measure');

% Cycle through each image and assign measurements
for i = 1:length(cells)
    nucleus = labels.Nucleus==cells(i);
    nucleus_erode = imerode(labels.Nucleus==cells(i),diskstrel(1));

    cytoplasm = (labels.Cell==cells(i)) &~nucleus;
    a = median(nfkb(nucleus));
    a_erode = median(nfkb(nucleus_erode));

    b = (mean(nfkb(cytoplasm))-distr2(1));
    if b<0
        b = 0;
    end
    b = b*numel(nfkb(cytoplasm));
    %disp(['nuc: ', num2str(a),'. cyto: ',num2str(b)])
    % Assign measurements
    CellMeasurements.NFkBdimNuclear(cells(i),ModuleData.iter) = a;
    CellMeasurements.NFkBdimNuclear_erode(cells(i),ModuleData.iter) = a_erode;
    CellMeasurements.NFkBdimCytoplasm(cells(i),ModuleData.iter) = b;
end

ModuleDataOut = ModuleData;

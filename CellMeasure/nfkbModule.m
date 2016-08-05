function [CellMeasurements, ModuleDataOut] = nfkbModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMODULE  measure nuclear fraction in auxiliary (fluorescent) image
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Grab NFkB image. If secondary AuxImage is defined, subtract it from NFkB image (AuxImage{1})

% Mode-balance 1st auxililiary image - bimodal distribution assumed
    if ~isfield(ModuleData,'distr')
        [AuxImages{1}, ModuleData.distr] = modebalance(AuxImages{1},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{1} = modebalance(AuxImages{1},2,ModuleData.BitDepth,'correct',ModuleData.distr);
    end


if ~isempty(AuxImages{2})
    % Mode-balance 2nd auxililiary image - bimodal distribution assumed
    if ~isfield(ModuleData,'distr2')
        [AuxImages{2}, ModuleData.distr2] = modebalance(AuxImages{2},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{2} = modebalance(AuxImages{2},2,ModuleData.BitDepth,'correct',ModuleData.distr2);
    end

    % Hard-coded correction factor (done by experiment)
    ind1 = regexp(parameters.ImagePath,'20..-..-');
    date1 = parameters.ImagePath(ind1:ind1+9);
    switch date1
        case '2012-11-05'
            weight = 0.07;
        case '2013-08-04'
            weight = 0.11;
        case '2013-08-08'
            weight = 0.11;        
        otherwise
            weight = 0.08;
    end
    % Correct image
    nfkb = AuxImages{1} - weight*AuxImages{2};
else
    nfkb = AuxImages{1};
end


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

% Cycle through each image and assign measurements
for i = 1:length(cells)
    bins = 0:4:(2^ModuleData.BitDepth)-1;
    nucleus = imerode(labels.Nucleus==cells(i),diskstrel(floor(parameters.MinNucleusRadius/5)));
    cytoplasm = (labels.Cell==cells(i)) &~nucleus;
    
    % Use median of nuclei compartment (mode would probably also work just fine)
    median_n = median(nfkb(nucleus));

    % Get secondary (higher) mode of cytoplasmic compartment
    n = hist(nfkb(cytoplasm),bins);
    thresh1 = quickthresh(nfkb,~cytoplasm,'none');
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

ModuleDataOut = ModuleData;

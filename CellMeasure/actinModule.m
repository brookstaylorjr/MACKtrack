

function [CellMeasurements, ModuleData] = actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ACTINMODULE identifies actin fibers in labeled cells, quantifying their orientation and scoring a cell by the 
% relativel level of organized fibers.
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%
%
% Methodology based on following sources:
% 1) Zemel A1, Rehfeldt F, Brown AE, Discher DE, Safran SA. Optimal matrix rigidity for stress fiber polarization in stem
% cells. Nat Phys. 2010 Jun 1;6(6):468-473.
% 2) Eltzner B, Wollnik C, Gottschlich C, Huckemann S, Rehfeldt F . The Filament Sensor for Near Real-Time Detection of
% Cytoskeletal Fiber Structures. PLOS ONE. 2015. 10(5): e0126346. https://doi.org/10.1371/journal.pone.0126346

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Misc. setup
measure_cc = label2cc(labels.Cell,0);
iteration  = ModuleData.iter;
actin_img = AuxImages{1};

% Parameters for fiber identification - assign by nuclear size range (correspond to roughly 10x, 20x, and 40x @ 2x2 bin)
gauss_std = 0.6*parameters.MinNucleusRadius; % Rough "best fit" made on a few test cases
steps = parameters.MinNucleusRadius+3;
min_size = round(parameters.MinNucleusRadius*2.5); % Fibers need > 1 nuclear diameter across


% Step 1: do eLoG processing, using rotated line Gaussians (from -pi/2 to pi/2)
g = gauss2D(gauss_std);
len = size(g,1);
f1 = [0 -1 0; -1 4 -1; 0 -1 0];
[x, y] = pol2cart(linspace(pi/2,-pi/2,steps), floor(len/2));
x = round(x(1:end-1)); y = round(y(1:end-1));
eLog = [];
origin = round(len/2);
for i = 1:(steps-1)
    r  = round(linspace(origin+y(i),origin-y(i),len));
    c  = round(linspace(origin+x(i),origin-x(i),len));
    rod = zeros(size(g));
    rod(sub2ind(size(g),r,c)) = g(sub2ind(size(g),r,c));
    eLog = cat(3,eLog,imfilter(imfilter(actin_img,rod,'symmetric'),f1,'symmetric'));    
end
max_line = max(eLog,[],3);
[~,idx] = sort(eLog,3,'descend');

% Step 2: mask the fiber image, then cleanup: preserve only objects with similar directions
mask = max_line>tsaithresh(max_line,labels.Cell==0);
angles = idx(:,:,1); 
angles(~mask) = nan;
mask1 = false(size(mask));
for i = 1:(steps)
    directions = circshift(1:steps,-(i-2),2);
    directions = directions(1:2);
    mask1 = mask1 | bwareaopen(ismember(angles,directions),min_size,8);
end
angles(~mask1) = nan;




% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'ActinScore')
    % Intensity-based measurement initialization
    CellMeasurements.ActinScore =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.ActinDist =  cell(parameters.TotalCells,parameters.TotalImages);
end



for i= 1: measure_cc.NumObjects
    vals = angles(measure_cc.PixelIdxList{i}); 
    vals(isnan(vals)) = [];
    n = histcounts(vals,0.5:(steps+0.5));
    
    % Score object; record score
    score = (max(n) - median(n))/sqrt(length(measure_cc.PixelIdxList{i}));
    CellMeasurements.ActinScore(i,iteration) = score;
    CellMeasurements.ActinDist{i,iteration} = n/sqrt(length(measure_cc.PixelIdxList{i}));
end



% (if a save directory is defined, e.g. by end_actinModule - i.e. don't output this in a live-cell context)
if isfield(ModuleData,'save_subdir')
    ModuleData.actin_dir = [save_subdir,filesep,'ActinModule',filesep];
    if ~exist(ModuleData.actin_dir,'dir');  mkdir(ModuleData.actin_dir); end
end


% Save a diagnostic image: show identified fiber angles overlaid on original
if isfield(ModuleData,'actin_dir')
    colormaps = loadcolormaps;
    tmp1 = actin_img;
    tmp1(tmp1==min(tmp1(:))) = []; tmp1(tmp1==max(tmp1(:))) = [];
    tmp1 = modebalance(tmp1,0, ModuleData.BitDepth,'display');               
    % Non-confluent case - set low saturation @ 3xS.D. below bg level
    if parameters.Confluence ~= 1
        pct = 90:.5:99;
        hi_val = prctile(tmp1,pct);       
        saturation_val = [-3 prctile(tmp1,1+findelbow(pct,hi_val))];
    else % Confluent case: unimodal distribution is foreground - use a different lower limit.
        saturation_val = [-4 prctile(tmp1(:),90)];
    end
    alpha = 0.55;
    baseImg = modebalance(actin_img,0,ModuleData.BitDepth);

    img = (baseImg - min(saturation_val))/ range(saturation_val);
    img(img<0) = 0 ; img(img>1) = 1;
    img = (img*255);
    base_cmap = colormaps.spectrum;
    cmap = [0 0 0 ; base_cmap(round(linspace(1,size(base_cmap,1),nanmax(angles(:)))),:)];
    img_overlay = angles;
    img_overlay(isnan(img_overlay))= 0;
    R = img*alpha + (1-alpha)*255*reshape(cmap(1+img_overlay(:),1),size(img)); 
    G = img*alpha + (1-alpha)*255*reshape(cmap(1+img_overlay(:),2),size(img)); 
    B = img*alpha + (1-alpha)*255*reshape(cmap(1+img_overlay(:),3),size(img)); 
    imwrite(uint8(cat(3,R,G,B)),[ModuleData.actin_dir,'ActinAngles-pos_',numseq(ModuleData.i,3),'.jpg'])
    
    
end

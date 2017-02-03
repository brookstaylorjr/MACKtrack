function handlesOut = loadImages(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% handlesOut = loadImages(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LOADIMAGES: use nucleus/cell image expressions to load up images for parameter 
% determination and testing
% 
% INPUTS:
% handles        master structure with object handles and parameters
%
% OUTPUTS:
% handlesOut     (updated) master structure with object handles and parameters
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% SETUP - clear all axes.
cla(handles.axes5A)
cla(handles.axes5B)
cla(handles.axes6A)
cla(handles.axes6B)
drawnow
parameters = handles.parameters;
locations = handles.locations;

% LOAD IMAGES
i = min(parameters.XYRange);
j = min(parameters.TimeRange);
nucleusFile = [locations.scope,parameters.ImagePath,eval(parameters.NucleusExpr)];

imfo = imfinfo([locations.scope,parameters.ImagePath,eval(parameters.NucleusExpr)]);
bit_depth = imfo.BitDepth;

nucOrig = checkread(nucleusFile,bit_depth);
parameters.ImageSize = [size(nucOrig,1),size(nucOrig,2)];
if ~strcmpi(parameters.ImageType,'none')
	cellFile = [locations.scope,parameters.ImagePath,eval(parameters.CellExpr)];
	cellOrig = checkread(cellFile,bit_depth);
    else
end


% NUCLEUS image
imfo = imfinfo(nucleusFile);
bit_depth = imfo.BitDepth;
baseImg = modebalance(nucOrig,1,bit_depth);
x = [-2 5];
baseImg(baseImg<x(1)) = x(1);
baseImg(baseImg>x(2)) = x(2);
baseImg = (baseImg - x(1))/(x(2)-x(1))*255;
nucDisp = uint8(round(baseImg));
imshow(nucDisp,'Parent',handles.axes5B)

% NUCLEUS OVERLAY (axes5B)
nucCircle1 = getnhood(diskstrel(parameters.MinNucleusRadius));
nucCircle2 = getnhood(diskstrel(parameters.MaxNucleusRadius));
diskNuc = padarray(nucCircle1,floor((parameters.ImageSize-size(nucCircle1))/2),0,'pre' );
diskNuc_temp = padarray(nucCircle2,floor((parameters.ImageSize-size(nucCircle2))/2),0,'pre' );
diskNuc = padarray(diskNuc,ceil((parameters.ImageSize-size(nucCircle1))/2),0,'post' ) + ...
    padarray(diskNuc_temp,ceil((parameters.ImageSize-size(nucCircle2))/2),0,'post' );
hold(handles.axes5B,'on')
handles.NucleusOverlay = imshow(label2rgb(diskNuc,'spring'),'Parent',handles.axes5B); 
set(handles.NucleusOverlay,'AlphaData',double(diskNuc>0)*0.6)
hold(handles.axes5B,'off')


% NUCLEUS EDGE PLOT: do initial noisecount search for parameter picking (axes5B)
nucSmooth = medfilt2(nucOrig,[parameters.MedianFilterSize, parameters.MedianFilterSize]);
nucEdgeHor = imfilter(nucSmooth,fspecial('sobel') /8,'replicate');
nucEdgeVert = imfilter(nucSmooth,fspecial('sobel')'/8,'replicate');
nucEdge = sqrt(nucEdgeHor.^2 + nucEdgeVert.^2);
levelStart = quickthresh(nucEdge,false(size(nucEdge)),'none');
imgSubsetNuc = imdilate(nucEdge>(levelStart), ones(80));
initialSearch = prctile(nucEdge(:),[5, 99]); % Use 5th and 99th percentile as starting point for nuclei
[noiseCount1, val1]  = noisecount(nucEdge,~imgSubsetNuc,initialSearch ,64);
set(handles.figure1,'CurrentAxes',handles.axes5A)
plot(handles.axes5A,val1,noiseCount1,'LineWidth',2,'Color',handles.orange)
set(handles.axes5A,'YTick',[],'FontSize',11)
% Keep starting value, unless it's way out of 5/95th prctile range.
nuc_edge = parameters.NucleusEdgeThreshold;
if (nuc_edge>(initialSearch(2)*2))||(nuc_edge<(initialSearch(1)/2))
    nuc_edge = mean(initialSearch);
    parameters.NucleusEdgeThreshold = nuc_edge;
    disp('Note: nuclear edge threshold was out of range; reset to middle of search range')
end
% Display line indicating initial search value
x_lims = [floor(min(initialSearch(1),nuc_edge)), ceil(max(initialSearch(2),nuc_edge))];
hold(handles.axes5A,'on')
handles.LineNuc1 = plot(handles.axes5A,ones(1,2)*nuc_edge,get(handles.axes5A,'YLim'),'Color',handles.blue);
hold(handles.axes5A,'off')
set(handles.slider5A,'Min',x_lims(1),'Max',x_lims(2),'Value',nuc_edge)
set(handles.edit5A,'String',num2str(nuc_edge))


% CELL OVERLAY (axes6B)
if ~strcmpi(parameters.ImageType,'none')
	% Phase contrast/DIC image: balance mode at 90.
	cellMin = cellOrig - min(cellOrig(:));
	[n, x] = hist(cellMin(:),2^bit_depth);
	imageMode = x(n==max(n));
	scale_1 = 90/imageMode(1);
	cellDisp = cellMin*scale_1;
	cellDisp(cellDisp>255)= 255;
	cellDisp = uint8(round(cellDisp));
	imshow(cellDisp,'Parent',handles.axes6B)
	% Cell overlays (cellCircle1/cellCircle2)
	radius1 = round(sqrt(parameters.MaxHoleSize/pi));
	cellCircle1 = getnhood(diskstrel(radius1));
	radius2 = round(sqrt(parameters.MinHoleSize/pi));
	cellCircle2 = getnhood(diskstrel(radius2));
	diskCell = padarray(cellCircle1,floor((parameters.ImageSize-size(cellCircle1))/2),0,'pre' );
	diskCell_temp = padarray(cellCircle2,floor((parameters.ImageSize-size(cellCircle2))/2),0,'pre' );
	diskCell = padarray(diskCell,ceil((parameters.ImageSize-size(cellCircle1))/2),0,'post' ) + ...
		padarray(diskCell_temp,ceil((parameters.ImageSize-size(cellCircle2))/2),0,'post' );
	hold(handles.axes6B,'on')
	handles.CellOverlay = imshow(label2rgb(diskCell,'winter'),'Parent',handles.axes6B); 
	set(handles.CellOverlay,'AlphaData',double(diskCell>0)*0.3)
	hold(handles.axes6B,'off')
    
    % CELL EDGE PLOT (axes6A)
    if ~strcmpi(parameters.ImageType,'fluorescence')

        horizontalEdge = imfilter(cellOrig,fspecial('sobel') /8,'replicate');
        verticalEdge = imfilter(cellOrig,fspecial('sobel')'/8,'replicate');
        edge_mag = sqrt(horizontalEdge.^2 + verticalEdge.^2);
        levelStart = quickthresh(edge_mag,false(size(edge_mag)),'none');
        imgSubsetCell = imdilate(edge_mag>(levelStart), ones(80));
        initialSearch = prctile(edge_mag(:),[5, 95]); % Use 5th and 95th percentile as starting point for cells
        [noiseCount2, val2]  = noisecount(edge_mag,~imgSubsetCell,initialSearch,64);
        plot(handles.axes6A,val2,noiseCount2,'LineWidth',2,'Color',[207 79 51]/255)
        set(handles.axes6A,'YTick',[],'FontSize',11)
        % Keep starting value, unless it's way out of 5/95th prctile range.
        cell_edge1 = min(parameters.CellSearchRange);
        if (cell_edge1<(initialSearch(1)/2))
            cell_edge1 = initialSearch(1);
            parameters.CellSearchRange = [initialSearch(1),max(parameters.CellSearchRange)];
            disp('Note: lower cell edge threshold was out of range; reset to bottom of search range')
        end
        cell_edge2 = max(parameters.CellSearchRange);
        if (cell_edge2>(initialSearch(2)*1.5))
            cell_edge2 = initialSearch(2);
            parameters.CellSearchRange = [min(parameters.CellSearchRange),initialSearch(2)];
            disp('Note: upper cell edge threshold was out of range; reset to top of search range')
        end
        % Plot lines indicating search range
        x_lims = [floor(min(initialSearch(1),cell_edge1)), ceil(max(initialSearch(2),cell_edge2))];
        hold(handles.axes6A,'on')
        handles.LineNoise1 = plot(handles.axes6A,ones(1,2)*cell_edge1,get(handles.axes6A,'YLim'),'Color',handles.blue);
        handles.LineNoise2 = plot(handles.axes6A,ones(1,2)*cell_edge2,get(handles.axes6A,'YLim'),'Color',handles.blue);
        hold(handles.axes6A,'off')
        set(handles.slider6A,'Min',x_lims(1),'Max',x_lims(2),'Value',cell_edge1)
        set(handles.slider6B,'Min',x_lims(1),'Max',x_lims(2),'Value',cell_edge2)
        set(handles.edit6A,'String',num2str(cell_edge1))
        set(handles.edit6B,'String',num2str(cell_edge2))
    end
end

% Save updated handles
handles.parameters = parameters;
handlesOut = handles;

% Temporary: throw loaded images into workspace
assignin('base','nucOrig',nucOrig)

% ========================================================================================

function [noiseCount, val] =noisecount(edges,dropPixels,searchRange,noiseSize)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% NOISECOUNT calculate noise pixels as a function of edge threshold
%
% edges         input edge magnitudes
% dropPixels    background area of images; normalizes over cell densities
% searchRange   min and max values that set threshold search range 
% noiseSize     maximum size of "noisy" pixel groups
%
% noiseCount    number of noise-sized object pixels at each value tested
% val           edge values tested
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% ________1.Take input edge image and scan over values, constructing noiseCount vector_____
val = min(searchRange): (max(searchRange)-min(searchRange))/100:max(searchRange);
edges(dropPixels) = 0;
noiseCount = zeros(size(val));

    for i = 1:length(val)
        mask = (abs(edges)>val(i));
        mask2 = bwareaopen(mask,noiseSize);
        sum1 = sum(mask(:)-mask2(:));
        noiseCount(i) = sum1;
    end

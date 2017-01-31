function [label_subobj, cc_obj] = edgesplit(label_in,split_img,split_img2,  p, mask_cut)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [label_subobj, label_obj] = edgesplit(label_in,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EDGESPLIT 
%
% SPLIT_IMG    either 1) Sobel edge magnitude (larger nuclei - use strong edges to split nuclei)
%              or 2) intensity watershed transform (smaller nuclei)
% SPLIT_IMG2    either 1) Inflection cut points (large nuclei)
%              or 2) original nucleus image (small nuclei)

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - SMALL NUCLEI: use strong edges to split regions - - - - - - - -
if p.MinNucleusRadius < 6
    cc_obj = label2cc(label_in,0);

    w1 = imdilate(split_img,ones(3));
    w1(label_in==0) = 0;
    pairs  = [w1(:),label_in(:)];
    [~,~,ic] = unique(pairs,'rows');
    label_subobj = reshape(ic,size(w1))-1;
    % Simplify objects in label1c to prevent bridgenuclei from hanging
    % [Count subobjects per larger object - cap @ 5]
    max_complexity = 6; % (Nuclei should already be mostly broken up by this point)
    get_obj = @(pxlist) unique(label_subobj(pxlist));
    obj_match = cellfun(get_obj, cc_obj.PixelIdxList,'UniformOutput',0)';
    complex_obj = find(cellfun(@length,obj_match)>max_complexity);
    % If complex objects are found, replace them with smoothed/re-watershedded image.
    n = 1;
    tmp_mask = imdilate(label_in>0,ones(3));
    while ~isempty(complex_obj)
        % Use larger smoothing kernel, recalculate watershed, and replace "complex" subregions as required
        nuc_smooth2 = imfilter(split_img2,gauss2D(0.33*n+p.MinNucleusRadius/4),'replicate'); % Gaussian filtered
        watershed2 = imdilate(watershedalt(nuc_smooth2, tmp_mask, 8),ones(3));
        for i = 1:length(complex_obj)
            subregion = cc_obj.PixelIdxList{complex_obj(i)};
            subregion_vals = unique(watershed2(subregion));
            lut = zeros(1,1+max(subregion_vals));
            lut(subregion_vals+1) = double(max(label_subobj(:)))+(1:length(subregion_vals));
            label_subobj(subregion) = lut(watershed2(subregion)+1);
        end
        % Re-count subobjects
        label_subobj = imclose(label_subobj,ones(2));
        get_obj = @(pxlist) unique(label_subobj(pxlist));
        obj_match = cellfun(get_obj, cc_obj.PixelIdxList,'UniformOutput',0)';
        complex_obj = find(cellfun(@length,obj_match)>max_complexity);
        n = n+1;
        label_subobj(label_in==0) = 0;
    end
    
% - - - - - - - LARGER NUCLEI: USE STRONG EDGES TO SUBDIVIDE REGIONS - - - - - - - 
else
    % Isolate high edge (>median) pixels per object
    tmp_cc = label2cc(label_in,0);
    get_high = @(pix) pix(split_img(pix) > prctile(split_img(pix),50));
    mask_vals = cell2mat(cellfun(get_high,tmp_cc.PixelIdxList,'UniformOutput',0));
    mask_hi = false(size(label_in));
    mask_hi(mask_vals) = 1;
    % Morphological cleanup and connection - dilate out, then shrink back result
    mask_hi = bwareaopen(mask_hi,round(p.MinNucleusRadius/2));
    mask_hi = ~bwareaopen(~mask_hi,round(p.MinNucleusRadius+1),4);
    BWconnect = bwmorph(mask_hi,'skel','Inf');
    for i = 3:2:round(sqrt(p.MinNucleusRadius))
        BWconnect = imdilate(BWconnect,ones(i));
        BWconnect = bwmorph(BWconnect,'skel','Inf');
        % Reduce back result of dilation
        for j = 1:floor(i/2)
            BWendpoints = bwmorph(bwmorph(BWconnect,'endpoints'),'shrink',Inf);
            BWconnect(BWendpoints) = 0;
        end
    end
    mask_hi = BWconnect;
    mask_obj = label_in>0;
    off_mask = mask_hi|BWconnect;
    % For "solid" objects (i.e. entirely high-edge), get borders of these directly and substitute
    solid_mask = imopen(mask_hi,diskstrel(p.NuclearSmooth));
    solid_borders = solid_mask&~imerode(solid_mask,ones(3));
    off_mask(solid_mask) = solid_borders(solid_mask);
    % Break up existing mask w/ newly-found borders -> use propagate subfcn to fill existing mask.
    mask_obj(off_mask) = 0;
    mask_obj = imopen(mask_obj,ones(2));
    min_area = floor(pi*(p.MinNucleusRadius)^2);
    mask_obj = bwareaopen(mask_obj,round(min_area/4),4);
    mask_all = label_in>0;
    label_subobj = IdentifySecPropagateSubfunction(double(bwlabel(mask_obj,4)),double(mask_all),mask_all,0.02);
    
    % Combine inflection point & strong edge data - see if a strong edge unabigiously connects two moderately-inflected pts.
    % (Filter out objects that are too small to be split further)
    filter_areas = @(pix) length(pix)<(2*min_area);
    nosmall = label2cc(label_in);
    nosmall.PixelIdxList(cellfun(filter_areas,nosmall.PixelIdxList)) = [];
    nosmall.NumObjects = length(nosmall.PixelIdxList);
    nosmall = labelmatrix(nosmall)>0;
    % Get borders from label1c, then identify if any lie a putative split point 
    tmp = label_subobj;
    tmp(~nosmall) = 0;
    tmp(tmp==0) = max(tmp(:))+1;
    border_mask = (tmp-imerode(tmp,ones(3)))>0;
    border_mask(label_subobj==0) = 0;
    endpt_val = bwmorph(border_mask,'endpoints').* split_img2;
    endpt_val(imdilate(mask_cut,ones(5))) = 0;
    endpt_val(endpt_val==0) = nan;
    branch_pts = bwmorph(bwmorph(border_mask,'skel',inf),'branchpoints');
    cc_branch = bwconncomp(border_mask,8);
    avg_inflect = @(pix) nanmean(endpt_val(pix));
    num_branch = @(pix) sum(branch_pts(pix));
    num_end = @(pix) sum(~isnan(endpt_val(pix)));
    check_avg = cellfun(avg_inflect,cc_branch.PixelIdxList);
    check_branch = cellfun(num_branch,cc_branch.PixelIdxList);
    check_ep = cellfun(num_end,cc_branch.PixelIdxList);
    
    mask_tmp = label_in>0;
    mask_tmp(cell2mat(cc_branch.PixelIdxList((check_avg>(p.NuclearInflection-30)) & (check_branch<1) & (check_ep==2))')) = 0;

    label_in2 = bwlabel(mask_tmp,4);
    label_subobj(label_in2==0) = 0;
    pairs  = [label_in2(:),label_subobj(:)];
    [~,~,ic] = unique(pairs,'rows');
    label_subobj = reshape(ic,size(label_in))-1;
    cc_obj = label2cc(label_in2);
    
end
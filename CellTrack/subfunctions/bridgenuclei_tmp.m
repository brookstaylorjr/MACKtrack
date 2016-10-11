label_in = diagnostics.nuclei.label1a;
label_in = imclose(label_in,ones(2));
cutoff.Area = [floor(pi*(p.MinNucleusRadius-1)^2) ceil(pi*(p.MaxNucleusRadius)^2)];
cutoff.Compactness = p.Compactness;
verbose = 0;
%%
tic
in_rprops = cell2mat(struct2cell(regionprops(label_in,'Area','Perimeter')))';
in_rprops = mat2cell(in_rprops,ones(1,size(in_rprops,1)),2);


% Combine connected objects and calculate regionprops
obj_cc = bwconncomp(label_in>0,4);
obj_rprops = cell2mat(struct2cell(regionprops(obj_cc,'Area','Perimeter')))';
obj_rprops = mat2cell(obj_rprops,ones(1,size(obj_rprops,1)),2);

% Assign initial objects (in label in) to groups of objects (in obj_cc)
get_obj = @(pxlist) unique(label_in(pxlist));
obj_match = cellfun(get_obj, obj_cc.PixelIdxList,'UniformOutput',0)';

% Get bounding boxes (and subsets) for all combined objects, to perform calculations on later
obj_bbox = struct2cell(regionprops(obj_cc,'BoundingBox'))';
get_subset = @(corners) label_in(max([1,floor(corners(2))]):min([size(label_in,1),ceil(corners(2)+corners(4))]),....
    max([1,floor(corners(1))]):min([size(label_in,2),ceil(corners(1)+corners(3))]));
label_subsets = cellfun(get_subset, obj_bbox,'UniformOutput',0);
drop_others = @(label, keep_idx) label .* uint16(ismember(label,keep_idx));
label_subsets = cellfun(drop_others, label_subsets, obj_match,'UniformOutput',0);

% Create checking functions
hard_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(1));
soft_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(2));

% PASS 1: Check intially-combined objects -> add in to pixels_out, drop from obj_cc
pass_1 = cellfun(hard_pass,obj_rprops);
pixels_out = obj_cc.PixelIdxList(pass_1);
if verbose
    for i = 1:length(pass_1)
        str_1 = ['obj containing [',num2str(obj_match{i}'),']. (a = ',num2str(obj_rprops{i}(1)),', c = ',...
                num2str((obj_rprops{i}(2).^2)/(4*pi*obj_rprops{i}(1))),')'];
        if pass_1(i)
            disp(['I. ADDED [',num2str(sum(pass_1(1:i))),']: ', str_1]);
        else
            disp(['FAILED ', str_1]);
        end
    end
end
remaining_obj = find(~pass_1);



% PASS 2: Of remaining (combined) objects, see if any is composed of subobjects that all hard_pass
pass_orig = cellfun(hard_pass, in_rprops);
all_pass = @(subobj) min(pass_orig(subobj));
pass_2 = cellfun(all_pass, obj_match(remaining_obj));
cc_in = label2cc(label_in,0);
if verbose
    add_list = cell2mat(obj_match(remaining_obj(pass_2)));
    for i = 1:length(add_list)
        disp(['II. ADDED [',num2str(i+length(pixels_out)),']: obj containing (all strict-passing) subobj [',...
            num2str(add_list(i)'),']'])
    end
end
pixels_out = [pixels_out(:); cc_in.PixelIdxList(cell2mat(obj_match(remaining_obj(pass_2))))];
remaining_obj = remaining_obj(~pass_2);


% PASS 3: soft pass on single/double objects
pass_3 = (cellfun(@length,obj_match(remaining_obj))<=2) & cellfun(soft_pass,obj_rprops(remaining_obj));
if verbose
    add_list = remaining_obj(pass_3);
    for i = 1:length(add_list)
        disp(['III. ADDED [',num2str(i+length(pixels_out)),']: obj containing [',num2str(obj_match{add_list(i)}),']. (a = ',...
            num2str(obj_rprops{add_list(i)}(1)),', c = ',...
            num2str((obj_rprops{add_list(i)}(2).^2)/(4*pi*obj_rprops{add_list(i)}(1))),')'])
    end
end
pixels_out = [pixels_out; reshape(obj_cc.PixelIdxList(remaining_obj(pass_3)),sum(pass_3),1)];
remaining_obj = remaining_obj(~pass_3);

% Drop out remaining single/double-subobject shapes
remaining_obj(cellfun(@length,obj_match(remaining_obj))==1) = [];


% Calculate neighbors for each original object (have to do 2 eroded versions to get b.g. as neighbor)
label_dil = imdilate(label_in,diskstrel(1));
tmp = label_in; tmp(tmp==0) = max(tmp(:))+1;
label_erd = imerode(tmp,diskstrel(1));
label_erd2 = imerode(label_in,diskstrel(1));

get_neighbors = @(set) unique([label_dil(set); label_erd(set); label_erd2(set)]);
neighbors = cellfun(get_neighbors,cc_in.PixelIdxList,'UniformOutput',0);


toc

%%

% REMAINING PASSES: for each of remaining component objects, 
for i = 1:50:length(remaining_obj);
    figure,imagesc(label_subsets{remaining_obj(i)})


end


toc

%% Calculate neighbors on original objects
tic
label_dil = imdilate(label_in,ones(3));
label_erd = imerode(label_in,ones(3));

get_neighbors = @(set) unique([label_dil(set); label_erd(set)]);
neighbors = cellfun(get_neighbors,cc_in.PixelIdxList,'UniformOutput',0);
toc


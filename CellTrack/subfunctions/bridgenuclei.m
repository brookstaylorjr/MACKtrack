function [label_out] = bridgenuclei(label_in,cutoff,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Check every object in a label mask, see if it can be bridged with anything else.
% Note: objects need to be touching each other in order to be bridged.
%
% old_label     label matrix to be merged
% cutoffs       morphological cutoffs: determine whether object is kept/merged/dropped. Cutoffs
%               must contain: area ([min max]), eccentricty (min), and compactness ([low high])
% verbose       boolean, 1 displays merging/testing output
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin <3
    verbose=0;
end

% Dilate label matrix over watershed borders
label_in = imclose(label_in,ones(2));

% Convert label_in to a conncomp struct. Calculate whether objects pass strict criteria (in isolation)
input_cc = label2cc(label_in,0);
tmp_props = regionprops(input_cc,'Area','Perimeter','MajorAxis','MinorAxis');
orig_obj = unique(label_in(label_in>0));
orig_pass = zeros(size(orig_obj));

for i = 1:length(orig_obj)
    % Test to see if fully merged object passes cutoffs 
    c = (tmp_props(orig_obj(i)).Perimeter.^2)/(4*pi*tmp_props(orig_obj(i)).Area);
    a = tmp_props(orig_obj(i)).Area;
    if (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && (c<cutoff.Compactness(1))
        orig_pass(i) = 1;
    end
end

if verbose
    var_idx = num2str(round(rand(1)*1000));
    assignin('base',['obj_pass_',var_idx],orig_pass)
end

% Combine connected objects and calculate regionprops
obj_cc = bwconncomp(label_in>0,4);
obj_rprops = regionprops(obj_cc,'Area','Perimeter');
label_out = zeros(size(label_in));
subobj_opened = cell(1);


% Iterate until every connected object from original has been evaluated
iter = 0; % make sure we don't get stuck.
while (obj_cc.NumObjects>0) && (iter<10000)
    obj_pass = 0;
    c = (obj_rprops(1).Perimeter.^2)/(4*pi*obj_rprops(1).Area);
    a = obj_rprops(1).Area;
    if verbose
        disp(['Testing obj containing [',num2str(unique(label_in(obj_cc.PixelIdxList{1}))'),']'])
        fprintf(['a = ',num2str(a),', c = ',num2str(c)])
    end
    % 1) STRICT check fully merged object
    if (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && (c<cutoff.Compactness(1))
        label_out(obj_cc.PixelIdxList{1}) = max(label_out(:))+1;
        obj_pass = 1;
        if verbose; disp(' -> added as new obj');  end
    else
        subobj = unique(label_in(obj_cc.PixelIdxList{1}));
        subobj(subobj==0) = [];
        % 2) STRICT check, all component objects
        if min(orig_pass(subobj))>0
            if verbose; disp('-> all component objects passed strict test; adding each as separate ohjects'); end
            for i = 1:length(subobj)
                label_out(input_cc.PixelIdxList{subobj(i)}) = max(label_out(:))+1;
            end
            obj_pass = 1;
            % 3) LOOSE check, full object  
        elseif (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && (c<cutoff.Compactness(2))
            label_out(obj_cc.PixelIdxList{1}) = max(label_out(:))+1;
            obj_pass = 1;
            if verbose; disp(' -> (passed looser criteria), added as new obj');  end
        end
    end  
        
        
    if ~obj_pass
        if verbose
            disp(' -> (fail)')
        end
        % FALSE; see if object contains >1 unique subobject
        if (length(subobj)>1)
            if (length(subobj)<12) % Keep object numbers down so combinatorial approach doesn't overwhelm                 
                % Step through all combinations of objects; exit if we find an object that fits critera
                flag_new = 0; % Reset flag
                for k = (length(subobj)-1):-1:1
                    cmb = nchoosek(subobj,k);
                    for row = 1:size(cmb,1)
                        subset_obj = cmb(row,:);
                        if verbose
                            disp(['Testing subobject(s) ',num2str(subset_obj)])
                        end
                        tmp_mask = false(size(label_in));
                        tmp_mask(cat(1,input_cc.PixelIdxList{subset_obj})) = 1;
                        tmp_cc = bwconncomp(tmp_mask,4);
                        % If combination yields connected obj, measure it and see if it passes cutoffs
                        if tmp_cc.NumObjects==1
                            tmp_rprops = regionprops(tmp_cc,'Area','Perimeter');
                            c = (tmp_rprops(1).Perimeter.^2)/(4*pi*tmp_rprops(1).Area);
                            a = tmp_rprops(1).Area;
                            if verbose
                                disp(['a = ',num2str(a),', c = ',num2str(c)])
                            end
                            % 1) STRICT test, fully merged object
                            obj_pass = 0;
                            if (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && (c<cutoff.Compactness(1))
                                label_out(tmp_mask) = max(label_out(:))+1;
                                obj_pass = 1;
                                if verbose; disp(['added as obj ',num2str(max(label_out(:))+1),'- adding subobjects back into stack']); end
                            % 2) STRICT check, all component objects
                            elseif min(orig_pass(subset_obj))>0
                                for m = 1:length(subset_obj)
                                    label_out(input_cc.PixelIdxList{subset_obj(m)}) = max(label_out(:))+1;
                                end
                                obj_pass = 1;
                                if verbose; disp(['all sub-objects passed; added them as obj ',num2str(max(label_out(:))-length(subset_obj)),' to ',num2str(max(label_out(:)))]); end
                            elseif (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && (c<cutoff.Compactness(1))
                                label_out(tmp_mask) = max(label_out(:))+1;
                                obj_pass = 1;
                                if verbose; disp(['(weaker criteria passed). Added as obj ',num2str(max(label_out(:))+1),'- adding subobjects back into stack']); end
   
                            end
                            % If object passed, take remaining subobjects in cluster and add (as new object(s)) back into stack
                            if obj_pass    
                                tmp_mask2 = false(size(label_in));
                                tmp_mask2(cat(1,input_cc.PixelIdxList{subobj})) = 1;
                                tmp_mask2(tmp_mask) = 0;
                                tmp_cc2 = bwconncomp(tmp_mask2,4);
                                tmp_rprops2 = regionprops(tmp_cc2,'Area','Perimeter');
                                obj_cc.PixelIdxList = cat(2,obj_cc.PixelIdxList,tmp_cc2.PixelIdxList);
                                obj_cc.NumObjects = obj_cc.NumObjects+tmp_cc2.NumObjects;
                                obj_rprops = [obj_rprops;tmp_rprops2];
                                flag_new = 1;
                                break
                            end
                        end
                    end
                    % drop out of k for loop if new object was made
                    if flag_new
                        break
                    end
                end % end k loop
            else % High # of subobjects in cluster: open the mask with a large element              
                flag_opened = 0; % if we haven't dealt with this cluster before, add new clusters back into stack
                for idx1 = 1:length(subobj_opened)
                    if isequal(subobj_opened{idx1},subobj)
                       flag_opened = 1; 
                    end
                end
                if ~flag_opened
                    tmp_mask2 = false(size(label_in));
                    tmp_mask2(cat(1,input_cc.PixelIdxList{subobj})) = 1;
                    radius1 = round(1.5*sqrt(cutoff.Area(1)/pi));
                    tmp_opened = imopen(tmp_mask2,diskstrel(radius1));
                    tmp_cc2 = bwconncomp(tmp_opened,4);
                    tmp_rprops2 = regionprops(tmp_cc2,'Area','Perimeter');
                    obj_cc.PixelIdxList = cat(2,obj_cc.PixelIdxList,tmp_cc2.PixelIdxList);
                    obj_cc.NumObjects = obj_cc.NumObjects+tmp_cc2.NumObjects;
                    obj_rprops = [obj_rprops;tmp_rprops2];
                    if verbose
                        disp(['imopen on subobjects [',num2str(subobj'),']'])
                    end
                    subobj_opened = cat(1,subobj_opened,subobj);
                end
            end
       end
    end
    % Done with object; delete it from structures
    obj_cc.PixelIdxList(1) = [];
    obj_cc.NumObjects = obj_cc.NumObjects-1;
    obj_rprops(1) = [];
    iter = iter+1;
    if iter==9999; warning('Max iterations reached'); end
end





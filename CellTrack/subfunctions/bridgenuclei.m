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

% dilate over watershed borders
label_in = imclose(label_in,ones(2));

% convert label_in to a conncomp struct
input_cc = label2cc(label_in);

% Combine connected objects and calculate regionprops
obj_cc = bwconncomp(label_in>0,4);
obj_rprops = regionprops(obj_cc,'Area','Perimeter','MajorAxis','MinorAxis');
label_out = zeros(size(label_in));
subobj_opened = cell(1);

% Iterate until everything's been checked in original
iter = 0; % make sure we don't get stuck.
while (obj_cc.NumObjects>0) && (iter<10000)
   % Test to see if fully merged object passes cutoffs 
    c = (obj_rprops(1).Perimeter.^2)/(4*pi*obj_rprops(1).Area);
    e = (obj_rprops(1).MajorAxisLength/obj_rprops(1).MinorAxisLength);
    a = obj_rprops(1).Area;
    if verbose
        disp(['Testing obj containing [',num2str(unique(label_in(obj_cc.PixelIdxList{1}))'),']'])
        fprintf(['e =',num2str(e), ', a = ',num2str(a),', c = ',num2str(c)])

    end
    
    if (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && ( (c<cutoff.Compactness(1)) ...
            || ((c<cutoff.Compactness(2))&&(e>cutoff.Eccentricity)) )
        % TRUE; add to label_out
        label_out(obj_cc.PixelIdxList{1}) = max(label_out(:))+1;
        if verbose
            disp(' -> added as new obj')
        end
    else
        if verbose
            disp(' -> (fail)')
        end
        % FALSE; see if object contains >1 unique subobject
        subobj = unique(label_in(obj_cc.PixelIdxList{1}));
        subobj(subobj==0) = [];
        if (length(subobj)>1)
            if (length(subobj)<12) % Keep object numbers down so combinatorial approach doesn't overwhelm                 
                % Step through all combinations of objects; exit if we find an object that fits critera
                flag_new = 0; % Reset flag
                for k = (length(subobj)-1):-1:1
                    cmb = nchoosek(subobj,k);
                    for row = 1:size(cmb,1)
                        if verbose
                            disp(['Testing subobject(s) ',num2str(cmb(row,:))])
                        end
                        tmp_mask = false(size(label_in));
                        tmp_mask(cat(1,input_cc.PixelIdxList{cmb(row,:)})) = 1;
                        tmp_cc = bwconncomp(tmp_mask,4);
                        % If combination yields connected obj, measure it and see if it passes cutoffs
                        if tmp_cc.NumObjects==1
                            tmp_rprops = regionprops(tmp_cc,'Area','Perimeter','MajorAxis','MinorAxis');
                            c = (tmp_rprops(1).Perimeter.^2)/(4*pi*tmp_rprops(1).Area);
                            e = (tmp_rprops(1).MajorAxisLength/tmp_rprops(1).MinorAxisLength);
                            a = tmp_rprops(1).Area;
                            if verbose
                                disp(['e =',num2str(e), ', a = ',num2str(a),', c = ',num2str(c)])
                            end
                            if (a < cutoff.Area(2)) && (a > cutoff.Area(1)) && ( (c<cutoff.Compactness(1)) ...
                                    || ((c<cutoff.Compactness(2))&&(e>cutoff.Eccentricity)) )
                                % If it passes, add object to label_out. Pull remaining objects and add them back to obj_cc
                                label_out(tmp_mask) = max(label_out(:))+1;
                                if verbose
                                    disp(['added as obj ',num2str(max(label_out(:))+1),'- adding other subobjects back into stack'])
                                end
                                tmp_mask2 = false(size(label_in));
                                tmp_mask2(cat(1,input_cc.PixelIdxList{subobj})) = 1;
                                tmp_mask2(tmp_mask) = 0;
                                tmp_cc2 = bwconncomp(tmp_mask2,4);
                                tmp_rprops2 = regionprops(tmp_cc2,'Area','Perimeter','MajorAxis','MinorAxis');
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
                    tmp_rprops2 = regionprops(tmp_cc2,'Area','Perimeter','MajorAxis','MinorAxis');
                    obj_cc.PixelIdxList = cat(2,obj_cc.PixelIdxList,tmp_cc2.PixelIdxList);
                    obj_cc.NumObjects = obj_cc.NumObjects+tmp_cc2.NumObjects;
                    obj_rprops = [obj_rprops;tmp_rprops2];
                    if verbose
                        disp(['imopen on subobects [',num2str(subobj'),']'])
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
end

function [SegmentData] = segmentSubset(well_subset, nuclear_channel, measurement_channels, measurement_type, image_dir, parameters, save_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [SegmentData] = segmentSubset(well_subset, measurement_channels, image_dir, output_dir, save_img)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SEGEMENTSUBSET performs primary object segmentation on all nuclear images in a specified folder (image_dir) that match 
% specified wells (well_subset, e.g. {'A01','A02'}). Mean intensity values in all measurement image (specified by
% measurement channels) will be measured in segmented nuclear objects, and output in SegmentData.
% 
% Filenames are constructed using default MicroIX conventions
% 
% INPUTS (required)
% well_subset           cell array of well names corresponding to cells to track
% nuclear_channel       string containing nuclear channel (e.g. 'w1')
% measurement_channels  cell array containing measurment channel(s) (e.g. {'w2','w3'})
% measurement_type      cell array containing measurement type, corresponding to channels (e.g {@mean, @mean})
% image_dir             directory where all original images are stored
% output_dir            directory where output files will be saved
% parameters            MACKtrack parameters structure
%
% INPUTS (optional)
% save_img              Boolean flag (1 or 0) that specifies whether an outlined image will be saved
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Default: don't save checking output
if nargin<7
    save_dir = '';
end


% Make output directory
if ~isempty(save_dir)
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
end

% Main loop: read images, segment, measure, and save output

% Concatenate full list of image names
image_names = struct2cell(dir(image_dir));
image_names = image_names(1,:)';

% Initialize output structures
SegmentData.image = {};
SegmentData.image_id = [];
SegmentData.cell_id = [];
SegmentData.measurements = [];
SegmentData.label_nuc = [];


for i = 1:length(well_subset)
    % Identify nuclear images using well names and nuclear channel
    nuc_id = ~cellfun(@isempty,strfind(image_names,['_',well_subset{i},'_']));
    nuc_id = nuc_id & ~cellfun(@isempty,strfind(image_names,['_',nuclear_channel]));
    % Drop anything with the name "thumb"
    nuc_id = nuc_id & cellfun(@isempty,strfind(image_names,'thumb'));
    nuc_id = find(nuc_id);
    if isempty(nuc_id)
        warning(['No corresponding images found for well ', well_subset{i}])
        continue
    end
        
    % Initialize measurements for this subset
    well_measurements = cell(size(nuc_id));
    well_image = cell(size(nuc_id));
    well_ids = cell(size(nuc_id));
    well_nuclei = cell(size(nuc_id));
    well_image_id = cell(size(nuc_id));
    
    if isempty(SegmentData.image_id)
        base_idx = 0;
    else
        base_idx = max(SegmentData.image_id);
    end
    
    % Cycle over all matching images - segment cells and measure.  
    for j = 1:length(nuc_id) 
        % Segment nuclear image
        tic
        nuc_orig = double(imread([image_dir,filesep,image_names{nuc_id(j)}]));
        output =  primaryID(nuc_orig,parameters,[]);
        output = nucleusID(nuc_orig,parameters,output);

        t = toc;
        % Save segmented nuclei (for diagnostic purposes)
        if exist(save_dir,'dir')
            tmp_img = nuc_orig;
            tmp_img = (tmp_img - prctile(tmp_img(:),3))/diff(prctile(tmp_img(:),[3 99]));
            tmp_img(tmp_img<0) = 0; tmp_img(tmp_img>1) = 1;
            tmp_img = tmp_img*255;
            tmp_mask = (imdilate(output.label_nuc,ones(3)) - output.label_nuc) > 0;
            R = tmp_img; R(tmp_mask) = 248*0.75 + R(tmp_mask)*0.25;
            G = tmp_img; G(tmp_mask) = 152*0.75 + G(tmp_mask)*0.25;
            B = tmp_img; B(tmp_mask) = 29*0.75 + B(tmp_mask)*0.25;
            tmp_img = uint8(cat(3,R,G,B));
            tmp_name = image_names{nuc_id(j)};
            imwrite(tmp_img,[save_dir,filesep,tmp_name(1:end-3),'.jpg'])
        end
        all_msg = ['Segmented ''', image_names{nuc_id(j)},'''. Elapsed time = ',num2str(t),' sec\n'];
        
        % (Make sure output mask is relabeled contiguously)
        nuc_cc = label2cc(output.label_nuc,'true');
        well_nuclei{j} = nuc_cc;
        % Measure objects in image and combine
        measurements = [];
        for k = 1:length(measurement_channels)
            % Get corresponding measurement image for each nuclear one
            measured_name = image_names{nuc_id(j)};            
            measured_name = [measured_name(1:strfind(measured_name,['_',nuclear_channel])),...
                measurement_channels{k}];
            measure_id = ~cellfun(@isempty,strfind(image_names,measured_name));
            measure_id = measure_id & cellfun(@isempty,strfind(image_names,'thumb')); % Drop anything with the name "thumb"
            measure_id = find(measure_id);
            if ~isempty(measure_id)
                all_msg = [all_msg, '(measurement col ',num2str(k),'): measuring ''', image_names{measure_id},' (',func2str(measurement_type{k}),')''\n'];
                measure_orig = double(imread([image_dir,filesep,image_names{measure_id}]));
            else
                measure_orig = nan(size(nuc_orig));
            end
            chan_measurement = zeros(nuc_cc.NumObjects,1);
            for m = 1:nuc_cc.NumObjects
                chan_measurement(m) = measurement_type{k}(measure_orig(nuc_cc.PixelIdxList{m}));
            end
            % Concatenate measurements together (as columns)
            measurements = cat(2,measurements,chan_measurement);
        end
        fprintf(all_msg)
        % Concatenate cells from grouped images together
        well_measurements{j} = measurements;
        well_ids{j} = (1:nuc_cc.NumObjects)';
        well_image{j} = image_names{nuc_id(j)};
        well_image_id{j} = repmat(j+base_idx,[size(measurements,1),1]);
    end
    
    % Concatenate in all the corresponding well measurements together
    SegmentData.image = cat(1,SegmentData.image,cellstr(cell2mat(well_image)));
    SegmentData.cell_id = cat(1,SegmentData.cell_id,cell2mat(well_ids));
    SegmentData.image_id = cat(1,SegmentData.image_id,cell2mat(well_image_id));
    SegmentData.measurements = cat(1,SegmentData.measurements,cell2mat(well_measurements));
    SegmentData.label_nuc = cat(1,SegmentData.label_nuc,well_nuclei);
end
% Add other supporting information
SegmentData.image_dir = image_dir;

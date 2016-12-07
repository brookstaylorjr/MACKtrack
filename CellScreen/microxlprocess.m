function Data = microxlprocess(image_dir, image_names, wells, save_subdir, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% measurements = microxlprocess(image_dir, image_names, wells, save_subdir, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MICROXLSPROCESS is called by screenLoop (in parallel) to parse an image directory, find
% wells/images corresponding condition, and segment/measure cells in those images.
%
% Micro XL imaging naming schema is assumed in identifying images - thumbnails are ignored.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Initialize substructures
Data.CellData = [];
Data.Measurements = struct;
Data.Images = {};
idx = 0;
for i = 1:length(wells)          
    % a) PER WELL: Identify all nuclear images in condition/well
    nuc_id = ~cellfun(@isempty,strfind(image_names,['_',wells{i},'_']));
    nuc_id = nuc_id & ~cellfun(@isempty,strfind(image_names,eval(parameters.NucleusMatch)));
    nuc_id = nuc_id & cellfun(@isempty,strfind(image_names,'thumb')); % Drop anything with the name "thumb"
    nuc_id = find(nuc_id);
    Data.Images = cat(1, Data.Images, image_names{nuc_id});
    if isempty(nuc_id)
        warning(['No corresponding images found for well ', wells{i}{k}])
        continue
    end 
    % b) PER IMAGE: Load, segment, save, and measure nuclear image
    for j = 1:length(nuc_id) 
        idx = idx+1; % Track total number of images used
        % 1) SEGMENT nuclear image
        tic
        nuc_orig = double(imread([image_dir,filesep,image_names{nuc_id(j)}]));
        output =  primaryID(nuc_orig,parameters,[]);
        output = nucleusID(nuc_orig,parameters,output);
        NuclearLabel = uint16(output.label_nuc);
        t1 = toc;
        % 2) SAVE diagnostic output and nuclear label matrix
        tic
        save_name = [save_subdir,filesep,'NuclearLabels',filesep,'NuclearLabel-',wells{i},'_',numseq(j,2),'.mat'];
        save(save_name,'NuclearLabel')
        save_name = [save_subdir,filesep,'SegmentedImages',filesep,'Segmentation-',wells{i},'_',numseq(j,2),'.jpg'];
        tmp_img = nuc_orig;
        tmp_img = (tmp_img - prctile(tmp_img(:),3))/diff(prctile(tmp_img(:),[3 99.2]));
        tmp_img(tmp_img<0) = 0; tmp_img(tmp_img>1) = 1;
        tmp_img = tmp_img*255;
        tmp_mask = (imdilate(output.label_nuc,ones(3)) - output.label_nuc) > 0;
        R = tmp_img; R(tmp_mask) = 17*0.75 + R(tmp_mask)*0.25;
        G = tmp_img; G(tmp_mask) = 255*0.75 + G(tmp_mask)*0.25;
        B = tmp_img; B(tmp_mask) = 58*0.75 + B(tmp_mask)*0.25;
        tmp_img = uint8(cat(3,R,G,B));
        imwrite(tmp_img,save_name)
        t2 = toc;
        % 3) MEASURE images according to selected modules
        tic
        CellMeasurements = struct;
        labels.Nucleus = NuclearLabel;
        labels.Cell = NuclearLabel;
        parameters.TotalCells = length(unique(NuclearLabel(NuclearLabel>0)));
        parameters.TotalImages = 1;
        ModuleData.BitDepth = parameters.BitDepth;
        ModuleData.iter = 1;
        if isfield(parameters,'Flatfield')
            ModuleData.Flatfield = parameters.Flatfield;
        end
        measure_names = '';
        
        for m = 1:length(parameters.ModuleNames)
            ModuleData.name = parameters.ModuleNames{m};
            AuxImages = cell(1,2);
            if  parameters.(ModuleData.name).Use == 1;                
                % Check/load/correct auxiliary images
                for aux = 1:2
                    % Check name
                    if aux==1; curr_expr = parameters.(ModuleData.name).ImageExpr;
                    else curr_expr = parameters.(ModuleData.name).(['ImageExpr',num2str(aux)]); end
                    try
                        measure_id = ~cellfun(@isempty,strfind(image_names,['_',wells{i},'_']));
                        measure_id = measure_id & ~cellfun(@isempty,strfind(image_names,eval(curr_expr)));
                        measure_id = measure_id & cellfun(@isempty,strfind(image_names,'thumb')); % Drop thumbs
                        measure_id = find(measure_id);
                        if length(measure_id)~=length(nuc_id)
                            disp('ERROR: nucleus image list: ')
                            disp(image_names(nuc_id))
                            disp('does not match up with measurement image list:')
                            disp(image_names(measure_id))
                            error(['Stopping screening in folder ', image_dir])
                        end
                        curr_name = [image_dir,filesep,image_names{measure_id(j)}];
                    catch
                        curr_name = '--';
                    end
                    % Check file and add it into AuxImages
                    if exist(curr_name,'file')
                        if ~isfield(ModuleData,'BitDepth')
                            ModuleData.BitDepth = imfo.BitDepth;
                        end
                        AuxImages{aux} = checkread(curr_name,ModuleData.BitDepth);
                        tmp_idx = regexp(curr_name,'[0-9]_w[0-9]');
                        tmp_idx = tmp_idx(1);
                        measure_names = [measure_names, ' + ' curr_name(tmp_idx-6:tmp_idx+3)];

                    end
                end
                % Call measurement function
                currentfn = str2func(ModuleData.name);
                [CellMeasurements, ModuleData] = ...
                    currentfn(CellMeasurements,parameters,labels, AuxImages, ModuleData);            
            end
        end
        % Combine measurements into a single master structure
        % - - - - - - CellData - - - - - - - -     
        % overall index | well index | cell index
        % - - - - - - - - - - - - - - - - - - - -
        Data.CellData = cat(1, Data.CellData,...
            [idx*ones(parameters.TotalCells,1),j*ones(parameters.TotalCells,1), (1:parameters.TotalCells)']);
        fnames = fieldnames(CellMeasurements);
        for n = 1:length(fnames)
            if ~isfield(Data.Measurements,fnames{n})
                Data.Measurements.(fnames{n}) = CellMeasurements.(fnames{n});
            else
                Data.Measurements.(fnames{n}) = cat(1,....
                    Data.Measurements.(fnames{n}),CellMeasurements.(fnames{n}));
            end
        end
        t3 = toc;

        % 4) DISPLAY status
        str = ['- - - - [',image_dir,'] - - - - - -'];
        str = sprintf([str,'\n', 'Segmentation (', image_names{nuc_id(j)}, ') - ',  num2str(t1),' sec ']);
        str = sprintf([str,'\n', 'Saving (', 'NuclearLabel-',wells{i},'_',numseq(j,2),'.mat', ') - ',  ...
            num2str(t2),' sec ']);
        str = sprintf([str,'\n', 'Measurement - ',  num2str(t3),' sec ']);
        str = sprintf([str,'\n (Measured: ',measure_names(3:end),')']);
        str = sprintf([str, '\n', '- - - - - - - - - -']);
        disp(str)  
    end         
end
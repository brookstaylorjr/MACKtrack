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
        warning(['No corresponding images found for well ', wells{i}])
        continue
    end
    %% Get corresponding list of cell images
    
    if strcmpi(parameters.ImageType,'fluorescence')
        if i==1
            cell_images = {};
        end        
        cell_id = ~cellfun(@isempty,strfind(image_names,['_',wells{i},'_']));
        cell_id = cell_id & ~cellfun(@isempty,strfind(image_names,eval(parameters.CellMatch)));
        cell_id = cell_id & cellfun(@isempty,strfind(image_names,'thumb')); % Drop anything with the name "thumb"
        cell_id = find(cell_id);
        cell_images = cat(1, cell_images, image_names{cell_id});
        % (quickdir returns images in alphabetical order, so well positions for nucleus/cell should be aligned)
    end
    %%
    
    % b) PER IMAGE: Load, segment, save, and measure nuclear image
    for j = 1:length(nuc_id) 
        idx = idx+1; % Track total number of images used
        % 1) SEGMENT nuclear image
        tic
        nuc_orig = double(imread([image_dir,filesep,image_names{nuc_id(j)}]));
        
        if strcmpi(parameters.ImageType,'fluorescence')
            cell_orig = double(imread([image_dir,filesep,image_names{cell_id(j)}]));
            output = fluorescenceID(cell_orig, parameters, nuc_orig);
            tmp = nucleusID(nuc_orig,parameters,output);
            output = combinestructures(tmp,output);
            output.nuclei = output.label_nuc; % (skip normal "check" function)
            tmp = fluorescenceSegment(output, cell_orig, parameters);
            output = combinestructures(tmp,output);
            NuclearLabel = uint16(output.nuclei);
            CellLabel = uint16(output.cells);
        else
            output =  primaryID(nuc_orig,parameters,[]);
            output = nucleusID(nuc_orig,parameters,output);
            NuclearLabel = uint16(output.label_nuc);
            CellLabel = NuclearLabel;
            cell_orig = nuc_orig;
        end
        t1 = toc;

        
        % 2) SAVE diagnostic output and nuclear (and cell, if defined) label matricies
        tic
        save_name = [save_subdir,filesep,'NuclearLabels',filesep,'NuclearLabel-',wells{i},'_',numseq(j,2),'.mat'];
        save(save_name,'NuclearLabel')
        if ~strcmpi(parameters.ImageType,'none')
            save_name = [save_subdir,filesep,'CellLabels',filesep,'CellLabel-',wells{i},'_',numseq(j,2),'.mat'];
            save(save_name,'CellLabel')
        end
        save_name = [save_subdir,filesep,'SegmentedImages',filesep,'Segmentation-',wells{i},'_',numseq(j,2),'.jpg'];

        tmp_img = cell_orig;
        tmp_img = (tmp_img - prctile(tmp_img(:),3))/diff(prctile(tmp_img(:),[3 99.2]));
        tmp_img(tmp_img<0) = 0; tmp_img(tmp_img>1) = 1;
        tmp_img = tmp_img*255;
        % Overlay nuclear borders as orange
        border1 = (imdilate(NuclearLabel,ones(3))-NuclearLabel)>0;
        R = tmp_img; R(border1) = R(border1)*0.25 + 0.75*248;
        G = tmp_img; G(border1) = G(border1)*0.25 + 0.75*152;
        B = tmp_img; B(border1) = B(border1)*0.25 + 0.75*29;
        if ~strcmpi(parameters.ImageType,'none')
            %  Overlay cell borders as light blue
            border1 = (imdilate(CellLabel,ones(3))-CellLabel)>0;
            R(border1) = R(border1)*0.15 + 0.85*118;
            G(border1) = G(border1)*0.15 + 0.85*180;
            B(border1) = B(border1)*0.15 + 0.85*203;
        end  
        tmp_img = uint8(cat(3,R,G,B));
        imwrite(tmp_img,save_name)
        t2 = toc;
        % 3) MEASURE images according to selected modules
        tic
        CellMeasurements = struct;
        labels.Nucleus = NuclearLabel;
        labels.Cell = CellLabel;
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
            [i*ones(parameters.TotalCells,1),idx*ones(parameters.TotalCells,1), (1:parameters.TotalCells)']);
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
        
        % Clean data: remove tmp measurmements
        if isfield(Data.Measurements,'tmp')
            Data.Measurements = rmfield(Data.Measurements,'tmp');  
        end

        % 4) DISPLAY status
        tmp1 = image_dir;
        tmp1(strfind(tmp1,'\')) = '/';
        str = ['- - - - [',tmp1,'] - - - - - -'];
        str = sprintf([str,'\n', 'Segmentation (', image_names{nuc_id(j)}, ') - ',  num2str(t1),' sec ']);     
        str = sprintf([str,'\n', 'Saving (', 'NuclearLabel-',wells{i},'_',numseq(j,2),'.mat', ') - ',  ...
            num2str(t2),' sec ']);
        str = sprintf([str,'\n', 'Measurement - ',  num2str(t3),' sec ']);
        str = sprintf([str,'\n (Measured: ',measure_names(3:end),')']);
        str = sprintf([str, '\n', '- - - - - - - - - -']);
        disp(str)  
    end         
end


% Add in cell image names (if applicable
if exist('cell_images','var')
    Data.Images = cat(2,Data.Images,cell_images);
end

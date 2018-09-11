 function Data = screenprocess(image_dir, image_names, wells, save_subdir, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% measurements = screenprocess(image_dir, image_names, wells, save_subdir, parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MICROXLSPROCESS is called by screenLoop (in parallel) to parse an image directory, find
% wells/images corresponding condition, and segment/measure cells in those images.
%
% See notes on supported naming schema in 'screenLoop' function
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Initialize substructures
Data.CellData = [];
Data.Measurements = struct;
Data.Images = {};
idx = 0;
for i = 1:length(wells)          
    % a) PER WELL: Identify all nuclear images in condition/well
    nuc_images = wellmatch(image_names,wells{i},eval(parameters.NucleusMatch),parameters.scope_type);
    Data.Images = cat(1, Data.Images, nuc_images);
    if isempty(nuc_images)
        warning(['No corresponding images found for well ', wells{i}])
        continue
    end
    %% Get corresponding list of cell images
    
    if strcmpi(parameters.ImageType,'fluorescence')
        if i==1
            all_cell_images = {};
        end        
        cell_images = wellmatch(image_names,wells{i},eval(parameters.CellMatch),parameters.scope_type);
        all_cell_images = cat(1,all_cell_images,cell_images);
        % (quickdir returns images in alphabetical order, so well positions for nucleus/cell should already be aligned)
    end
    %%
    
    % b) PER IMAGE: Load, segment, save, and measure nuclear image
    for j = 1:length(nuc_images) 
        idx = idx+1; % Track total number of images used
        
        % 1) SEGMENT nuclear image (and cell image, if applicable)
        tic
        nuc_orig = double(imread([image_dir,filesep,nuc_images{j}]));
        
        if strcmpi(parameters.ImageType,'fluorescence')
            cell_orig = double(imread([image_dir,filesep,cell_images{j}]));
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
        R = tmp_img; R(border1) = R(border1)*0.25 + 0.75*0;
        G = tmp_img; G(border1) = G(border1)*0.25 + 0.75*195;
        B = tmp_img; B(border1) = B(border1)*0.25 + 0.75*169;
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
        ModuleData.save_subdir = save_subdir; % in case it's required for additional diagnostic output
        ModuleData.i = [wells{i},'_',numseq(j,2)];
        ModuleData.BitDepth = parameters.BitDepth;
        ModuleData.iter = 1;
        if isfield(parameters,'Flatfield')
            ModuleData.Flatfield = parameters.Flatfield;
        end
        measure_names = '';
        
        for m = 1:length(parameters.ModuleNames)
            ModuleData.name = parameters.ModuleNames{m};
            AuxImages = cell(1,3);
            if  parameters.(ModuleData.name).Use == 1;                
                % Check/load/correct auxiliary images
                for aux = 1:3
                    % Check name
                    if aux==1; curr_measure = parameters.(ModuleData.name).ImageExpr;
                    else curr_measure = parameters.(ModuleData.name).(['ImageExpr',num2str(aux)]); end
                    try                    
                        measure_images = wellmatch(image_names,wells{i},eval(curr_measure),parameters.scope_type);
                        if ~isempty(measure_images) && (length(measure_images)~=length(nuc_images))
                            disp(['NOTE: number of nuclear images for well ', wells{i}, '(',num2str(length(nuc_images)),...
                                ' sites found) does not match up with # of measurement images (',...
                                num2str(length(measure_images)),' sites) - may cause measurement error!'])
                        end
                        curr_name = [image_dir,filesep,measure_images{j}];
                    catch
                        curr_name = '--';
                    end
                    % Check file and add it into AuxImages
                    if exist(curr_name,'file')
                        if ~isfield(ModuleData,'BitDepth')
                            ModuleData.BitDepth = imfo.BitDepth;
                        end
                        
                        AuxImages{aux} = checkread(curr_name,ModuleData.BitDepth,1,1);
                        if ~isequal(size(AuxImages{aux}), size(labels.Nucleus))
                            AuxImages{aux} = zeros(size(labels.Nucleus));
                            warning([curr_name, 'is invalid - replacing with blank image'])
                        end
                        measure_names = [measure_names, ' + ' , eval(curr_measure)];

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
        tmp1 = namecheck(image_dir);
        tmp1(strfind(tmp1,'\')) = '/'; % STUPID WINDOWS BACK SLASHES!
        str = ['- - - - [',tmp1,'] - - - - - -'];
        str = sprintf([str,'\n', 'Segmentation (', nuc_images{j}, ') - ',  num2str(t1),' sec ']);     
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
    Data.Images = cat(2,Data.Images,all_cell_images);
end

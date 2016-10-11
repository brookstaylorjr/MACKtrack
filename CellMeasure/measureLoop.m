function [] = measureLoop(xy, parameters) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = measureLoop(xy, parameters) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MEASURELOOP is called in a parfor loop inside MACKmeasure.mat - it measures the full
% sequence of images from on XY position.
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Form xy path, initialize structures for measurements and module information              
parameters.XYDir = [parameters.locations.data,filesep,parameters.SaveDirectory,filesep,'xy',num2str(xy),filesep];   
CellMeasurements = struct;      
ModuleData = struct;

% Load and add CellData field to CellMeasurements
if exist([parameters.XYDir,'CellData.mat'],'file')
    load([parameters.XYDir,'CellData.mat']);
    parameters.TotalCells = length(CellData.FrameIn);
    i = xy;

     % Form CellMeasurements.CellData (for filtering and display)
     CellMeasurements.CellData = ...
         [i*ones(parameters.TotalCells,1),... % 1) xy position of cell
         (1:parameters.TotalCells)',... % 2) Cell ID number
         CellData.FrameIn,... % 3) First appearance of cell
         CellData.FrameOut,... % 4) Last appearance of cell
         CellData.Parent,... % 5) ID number of parent
         CellData.Edge];   % 6) Boolean of whether cell touches image edge during lifetime


    % Inner loop: cycle time points (in each xy)
    tic
    for iter = 1:parameters.TotalImages
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Run measurement modules, adding new fields to CellMeasurements
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load label matricies ('CellLabel' and 'NuclearLabel')
        if exist([parameters.XYDir,'CellLabels'],'dir')
            load([parameters.XYDir,'CellLabels',filesep,'CellLabel-',numseq(iter,4),'.mat'])
            labels.Cell = CellLabel;
        end
        load([parameters.XYDir,'NuclearLabels',filesep,'NuclearLabel-',numseq(iter,4),'.mat'])
        labels.Nucleus = NuclearLabel; 
        % Make sure that there aren't extra cells or nuclei around
        if isfield(labels,'Cell')
            cells = unique(labels.Cell(:));
            nuclei = unique(labels.Nucleus(:));
            cells(cells==0) = [];
            nuclei(nuclei==0) = [];
            nocell = nuclei(~ismember(nuclei,cells));
            nonuc = cells(~ismember(cells,nuclei));
            labels.Nucleus(ismember(labels.Nucleus,nocell)) = 0;
            labels.Cell(ismember(labels.Cell,nonuc)) = 0;

            if (~isempty(nocell))||(~isempty(nonuc))
                disp(['Nucleus/cell(s) [', num2str([nonuc', nocell']),'] partially missing', ...
                    ' (xy ',num2str(i),' frame ',num2str(iter),')'])
            end
        end
        
        j = parameters.TimeRange(iter);
        

        % Cycle through module names- construct name, load aux image, and call
        for m = 1:length(parameters.ModuleNames)
            ModuleData.name = parameters.ModuleNames{m};
            ModuleData.iter = iter;
            ModuleData.i = i;
            ModuleData.j = j;
            AuxImages = cell(1,2);
            if  parameters.(ModuleData.name).Use == 1;                
                % Check/load/correct auxiliary images
                for aux = 1:2
                    % Check name
                    if aux==1
                        curr_expr = parameters.(ModuleData.name).ImageExpr;
                    else
                        curr_expr = parameters.(ModuleData.name).(['ImageExpr',num2str(aux)]);
                    end
                    try
                        curr_name = [parameters.locations.scope, filesep,parameters.ImagePath, filesep, eval(curr_expr)];
                    catch
                        curr_name = '--';
                    end
                    % Check file and add it into AuxImages
                    if exist(curr_name,'file')
                        if ~isfield(ModuleData,'BitDepth')
                            imfo = imfinfo(curr_name);
                            ModuleData.BitDepth = imfo.BitDepth;
                        end
                        AuxImages{aux} = checkread(curr_name,ModuleData.BitDepth);   
                    end
                end
                
                % Call measurement function
                currentfn = str2func(ModuleData.name);
                [CellMeasurements, ModuleData] = currentfn(CellMeasurements,parameters,labels, AuxImages, ModuleData);            
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display status on every 10th run
        if (mod(iter,10)==0)
            t_ten = toc;
            disp(['Measured xy ',num2str(xy),', frames ', num2str(iter-9), ' to ', num2str(iter),': ', num2str(t_ten),' sec'])
            tic;
        end
    end

    % Delete temporary fields, concatenate data from each XY position
    if isfield(CellMeasurements,'tmp')
        CellMeasurements = rmfield(CellMeasurements,'tmp');  
    end
    
    % Save CellMeasurements
    save([parameters.XYDir,'CellMeasurements.mat'], 'CellMeasurements')
else
    disp(['Skipping XY ',num2str(xy)])
end
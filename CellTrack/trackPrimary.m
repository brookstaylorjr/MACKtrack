function [] = trackPrimary(parameters,xyPos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [] = trackPrimary(parameters,xyPos)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% TRACKPRIMARY takes a time-lapse image series and converts into segmented, tracked data. 
% Output cellular trajetories are saved as successive label matricies. It is identical to
% TRACKLOOP, but finds only nuclei (primary objects)
%
% Main subfunctions
% phaseID.m/dicID.m, nucleusID.m, trackNuclei.m. dicCheck.m, phaseSegment.m/dicSegment.m
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
parameters.debug = 0;

% SETUP: define options, initialize structures for images/masks/label matricies
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')

images = struct;
tocs = struct;

 % Get image bit depth
i = xyPos;
j =  parameters.TimeRange(1);
imfo = imfinfo([locations.scope,parameters.ImagePath,eval(parameters.NucleusExpr)]);
bit_depth = imfo.BitDepth;

% Make save directories/image stacks
outputDirectory = [locations.data,filesep, parameters.SaveDirectory,filesep,'xy',num2str(xyPos),filesep];
mkdir(outputDirectory)
mkdir([outputDirectory,'NuclearLabels'])    
mkdir([outputDirectory,'SegmentedImages'])

% Check to make sure time vector is long enough
if length(parameters.TimeRange) < parameters.StackSize
   error('Time vector is too short for stack size, aborting.')
end

% Save tracking/memory checking text output
fid = fopen([outputDirectory,'decisions.txt'],'w','n','UTF-8');
fwrite(fid, sprintf(['Tracking/checking decisions for xy pos ',num2str(xyPos),':\n']));
fclose(fid);


% Loop all time points
for cycle = 1:length(parameters.TimeRange)
    tic
    trackstring = '';
    i = xyPos;
    j =  parameters.TimeRange(cycle);
       
    % NUCLEAR IDENTIFICATION
    tic
    nucName1 = eval(parameters.NucleusExpr);
    images.nuc = checkread([locations.scope,parameters.ImagePath,nucName1],bit_depth,1,parameters.debug);
    data = primaryID(images.nuc,parameters);
    tocs.NucMasking = toc;
    
    % Update stacks/structs with each iteration      
    % After fill loops, empty bottom on update
    if cycle>parameters.StackSize
        future(1) = [];
    end
    % Concatenate new information into 'future' queue
    if cycle==1
        future = data;
    else
        future = cat(1,future,data);
    end
   
   if cycle >= parameters.StackSize
       % Bookkeeping (indicies), initialization for tracking 
        saveCycle = cycle-parameters.StackSize+1; % Value assigned to CellData and tracked label matricies
        j = parameters.TimeRange(saveCycle); % Number of the input image corresponding to the BOTTOM of stack
        % Read 1st auxiliary image - corresponding to bottom of the stack
        images.bottom = checkread([locations.scope,parameters.ImagePath,eval(parameters.CellExpr)],...
            bit_depth,1,parameters.debug);       
        
        % TRACKING: Initialize CellData (blocks and CellData) when queue is full, then track nuclei
        tic
        if cycle == parameters.StackSize
            [CellData, future] = initializeCellData(future,parameters);
        else
            trackstring = [trackstring,'\n- - - Cycle ',num2str(saveCycle),' - - -\n'];
            [tmpstring, CellData, future] =  evalc('trackNuclei(future, CellData, saveCycle, parameters)');
            trackstring = [trackstring,tmpstring];
        end
        tocs.Tracking = toc;    
        past = future(1);
        % Final CellData cleanup: mark cells that touch the edge of image - - - -
        edgeCheck = past(1).nuclei;
        edgeCells = unique([edgeCheck(1,:),edgeCheck(end,:),edgeCheck(:,1)',edgeCheck(:,end)']);
        edgeCells(edgeCells==0) = [];
        CellData.Edge(edgeCells) = 1;
        
        % SAVING (label mats, segmentated images, decisions.txt)
        tic
        % Save nuclear labels
        NuclearLabel = uint16(past(1).nuclei);
        save([outputDirectory,'NuclearLabels',filesep,'NuclearLabel-',numseq(saveCycle,4),'.mat'], 'NuclearLabel')
        % Save composite 'Segmentation' image
        alpha = 0.30;
        saveFig(images.bottom,NuclearLabel,NuclearLabel,[],bit_depth,[outputDirectory,'SegmentedImages',filesep,'Segmentation-',numseq(saveCycle,4),'.jpg'],  alpha)

        tocs.Saving = toc;
        % Save decisions.txt
        fid = fopen([outputDirectory,'decisions.txt'],'a','n','UTF-8');
        fwrite(fid, sprintf(trackstring));
        fclose(fid);
   end
    
    
    
    % - - - - Display progress/times taken for mask/track/segment - - - -
    name1 = parameters.SaveDirectory;
    if strcmp(name1(end),filesep)
        name1 = name1(1:end-1);
    end
    seps = strfind(name1,filesep);
    seps = seps(end);
    name1 = name1(seps+1:end);
    str = '\n- - - - - - - - - - - - - - -';
    if cycle < parameters.StackSize 
        str = sprintf([str, '\n', name1, ' - XY ', num2str(xyPos),', Fill Cycle ', num2str(cycle)]);
    else
        str = sprintf([str, '\n', name1, ' - XY ', num2str(xyPos),', Save Cycle ', num2str(saveCycle)]);
    end
        str = sprintf([str, '\n', 'Nucleus Image: ',nucName1]);
    n = fieldnames(tocs);
    for k = 1:length(n)
        str = sprintf([str '\n', n{k},'- ',num2str(tocs.(n{k})),' sec']);
    end
    fprintf([str, '\n'])

end

% Save CellData
save([outputDirectory,'CellData.mat'],'CellData')


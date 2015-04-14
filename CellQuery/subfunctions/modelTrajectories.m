function handlesOut = modelTrajectories(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MODELTRAJECTORIES  performs statistical modeling/clustering on measurement data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Columns in CellData: 
% 1) xyPosition
% 2) index in xy pos
% 3) FrameIn 
% 4) FrameOut 
% 5) Parent 
% 6) EdgeCell

% Clear option-specific fields 
if isfield(handles.Export,'Gaussian')
    handles.Export = rmfield(handles.Export,'Gaussian');
end

if isfield(handles.Export,'Linkage')
    handles.Export = rmfield(handles.Export,'Linkage');
end


% Switch by grouping/modeling type
switch handles.Options.Grouping
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case get(handles.radiobutton8A,'String') % No sorting
        handles.Export.GroupingVector = (1:size(handles.Export.Measurement1,1))';
        sortorder = 1:size(handles.Export.Measurement1,1); % No clustering, sortorder is just 1:n cells
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case get(handles.radiobutton8B,'String') % By lineage
        % Lineage case: cycle through cells, grouping children with parents
        handles.Export.GroupingVector = zeros(size(handles.Export.Measurement1,1),1);
        for i = 1:length(handles.Export.GroupingVector)
            if handles.Export.CellData(i,5)==0
                handles.Export.GroupingVector(i) = max(handles.Export.GroupingVector)+1;
            else
                % Copy parent's information into child, then assign child to parent's group
                parentInd = (handles.Export.CellData(:,1)==handles.Export.CellData(i,1)) & (handles.Export.CellData(:,2) == handles.Export.CellData(i,5));
                if sum(parentInd)==1
                    handles.Export.GroupingVector(i) = handles.Export.GroupingVector(parentInd);
                    handles.Export.Measurement1(i,1:handles.Export.CellData(parentInd,4)) =... % Copy parent's Measurement1 information
                        handles.Export.Measurement1(parentInd,1:handles.Export.CellData(parentInd,4));
                    handles.Export.Measurement2(i,1:handles.Export.CellData(parentInd,4)) =... % Copy parent's Measurement2 information
                        handles.Export.Measurement2(parentInd,1:handles.Export.CellData(parentInd,4));                    
                else
                    handles.Export.GroupingVector(i) = max(handles.Export.GroupingVector)+1;              
                end
            end
        
        end
        [handles.Export.GroupingVector,sortorder] = sort(handles.Export.GroupingVector,'ascend');
   
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case get(handles.radiobutton8C,'String') % Gaussian Mix (1-D)
        % Check Measure2Grouping flag, select measurement to model/group
        if handles.Options.Measure2Grouping
            measurement = handles.Export.Measurement2;
        else
            measurement = handles.Export.Measurement1;
        end
        
        handles.Export.GroupingVector = zeros(size(measurement,1),1);
        Gaussian.Cluster = nan(size(handles.Export.Measurement1));
        warning off all
        disp('1-D Gaussian modeling/clustering...')
        % Make Gaussian distribution/assign cells to cluster at each separate timepoint
        for i = 1:size(measurement,2)
            Gaussian.Dist{i} = gmdistribution.fit(measurement(:,i),handles.Options.Subpopulations); % Use user-set number of subpopulations
            % Order Gaussian mixture components by mean and cluster
            [~,sortind] = sort(Gaussian.Dist{i}.mu,'ascend');
            % Need to go from sorted index to ranking
            [~, sortorder] = sort(sortind,'ascend');
            Gaussian.Cluster(:,i) = cluster(Gaussian.Dist{i},measurement(:,i));
            tmpcol = Gaussian.Cluster(:,i);
            tmpcol(~isnan(tmpcol)) = sortorder(tmpcol(~isnan(tmpcol)));
            Gaussian.Cluster(:,i) = tmpcol;
            if mod(i,25) == 0
                disp(['t = ',num2str(i)]);
            end
            
        end

        disp('Rearranging data for display...')
        % GroupingVector: sort by most common cluster, then length
        for i = 1:length(handles.Export.GroupingVector)
            row_tmp = Gaussian.Cluster(i,:);
            handles.Export.GroupingVector(i) = mean(row_tmp(~isnan(row_tmp)));
        end
        [handles.Export.GroupingVector,sortorder] = sort(handles.Export.GroupingVector);
        handles.Export.GroupingVector = round(handles.Export.GroupingVector);
        
        for i = 1:size(Gaussian.Cluster,2)
            Gaussian.Cluster(:,i) = Gaussian.Cluster(sortorder,i);
        end
        
        % Save Gaussian model info (for hist overlay)
        handles.Export.Gaussian = Gaussian;
        warning on all % Turn warnings back on
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case get(handles.radiobutton8D,'String') % Hierarchial (Euclidean)
        % Check Measure2Grouping flag, select measurement to model/group
        if handles.Options.Measure2Grouping
            measurement = handles.Export.Measurement2;
        else
            measurement = handles.Export.Measurement1;
        end
     
        % Compute point-by-point "median trajectory" and calculate differences
        median_traj = zeros(1,size(measurement,2));
        for i = 1:size(measurement,2)
            tmpcol = measurement(:,i);
            median_traj(i) = median(tmpcol(~isnan(tmpcol)));
            tmpcol(isnan(tmpcol)) = median_traj(i);
            measurement(:,i) = tmpcol;
        end
        median_diff = zeros(size(measurement));
        for i = 1:size(measurement,1)
            median_diff(i,:) = (measurement(i,:) - median_traj);

        end
        median_diff(isnan(median_diff))= 0;
        % Do clustering/ordering based on Euclidean distance 
        p  = pdist(median_diff,'euclidean');
        handles.Export.Linkage = linkage(p); 
        
        % Get sortorder from dendrogram (Supress initial dendrogram output)
        tmpfig = figure('Visible','off');
        figure(tmpfig),[~,~,perm] = dendrogram(handles.Export.Linkage,0);
        close(tmpfig)
        sortorder = perm(end:-1:1);
                
        % Could assign to clusters, but clustering gets dominated by outlier cells
        handles.Export.GroupingVector = (1:size(handles.Export.Measurement1,1))';
        %handles.Export.GroupingVector = cluster(Z,'maxclust',handles.Options.Subpopulations)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case get(handles.radiobutton8E,'String') % Gaussian (2-D)
        % Pull both measurements
        measurement1 = handles.Export.Measurement1;
        measurement2 = handles.Export.Measurement2;
        warning off all
        
        handles.Export.GroupingVector = zeros(size(measurement1,1),1);
        Gaussian.Cluster = nan(size(handles.Export.Measurement1));
        disp('2-D Gaussian modeling/clustering...')
        % Make Gaussian distribution/assign cells to cluster at each separate timepoint
        for i = 1:size(measurement1,2)
            Gaussian.Dist{i} = gmdistribution.fit([measurement1(:,i),measurement2(:,i)],handles.Options.Subpopulations); % Use user-set number of subpopulations
            % Order Gaussian mixture components by mean and cluster
            [~,sortind] = sort(Gaussian.Dist{i}.mu(:,1),'ascend');
            % Need to go from sorted index to ranking
            [~, sortorder] = sort(sortind,'ascend');
            Gaussian.Cluster(:,i) = cluster(Gaussian.Dist{i},[measurement1(:,i),measurement2(:,i)]);
            tmpcol = Gaussian.Cluster(:,i);
            tmpcol(~isnan(tmpcol)) = sortorder(tmpcol(~isnan(tmpcol)));
            Gaussian.Cluster(:,i) = tmpcol;
            if mod(i,25) == 0
                disp(['t = ',num2str(i)]);
            end
            
        end

        disp('Rearranging data for display...')
        % GroupingVector: sort by most common cluster, then length
        for i = 1:length(handles.Export.GroupingVector)
            row_tmp = Gaussian.Cluster(i,:);
            handles.Export.GroupingVector(i) = mean(row_tmp(~isnan(row_tmp)));
        end
        [handles.Export.GroupingVector,sortorder] = sort(handles.Export.GroupingVector);
        handles.Export.GroupingVector = round(handles.Export.GroupingVector);

        for i = 1:size(Gaussian.Cluster,2)
            Gaussian.Cluster(:,i) = Gaussian.Cluster(sortorder,i);
        end
        
        % Save Gaussian model info (for hist/scatterplot overlay)
        handles.Export.Gaussian = Gaussian;
        warning on all        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end

% Save sort order and assign handlesOut
handles.Export.Order = sortorder;
handlesOut = handles;
        
        
        
        
        
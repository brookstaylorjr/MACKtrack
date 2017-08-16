function [by_condition, by_well, by_image] = restructuredata(AllData, measure_name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [by_condition, by_well, by_image] = structuredata(AllData, measure_name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% RESTRUCTUREDATA splits a single-cell measurement (e.g. 'MeanNuc1') into 3 different cell matricies, organized:
% 1) condition
% 2) WELLS within a given condition
% 3) IMAGES within a given condition
%
% Example: [xdata_by_condition, xdata_by_well, xdata_by_image] = restructuredata(AllData,'MeanNuc2');
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Collect data together
all_cond = fieldnames(AllData);
by_condition = cell(size(all_cond));
by_well = cell(size(all_cond));
by_image = cell(size(all_cond));

for i = 1:length(all_cond)
    if ~isfield(AllData.(all_cond{i}).Measurements,measure_name)
        error(['Error: Measurement field ''', measure_name,''' does not exist in one or more conditions.'])
    end
    by_condition{i} =  AllData.(all_cond{i}).Measurements.(measure_name)(:);

    % Group by well
    wells = unique(AllData.(all_cond{i}).CellData(:,1));
    by_well{i} = cell(size(wells));
    for j = 1:length(wells)
        by_well{i}{j} = AllData.(all_cond{i}).Measurements.(measure_name)(AllData.(all_cond{i}).CellData(:,1)==wells(j));
    end
    
    % Group by image
    images = unique(AllData.(all_cond{i}).CellData(:,2));
    by_image{i} = cell(size(images));
    for j = 1:length(images)
        by_image{i}{j} = AllData.(all_cond{i}).Measurements.(measure_name)(AllData.(all_cond{i}).CellData(:,2)==images(j));
    end

end
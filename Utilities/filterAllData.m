function AllData_out = filterAllData(AllData, filter_field, filter_function)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% AllData_out = filterAllData(AllData_in, filter_function)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FILTERALLDATA applies a supplied filtering function (e.g. @(x) x>0) to an input filter field (string, e.g.
% 'MeanNuc1'). Any cells that meet this criteria (in above example, with negative values) will be filtered out.
%
%
% INPUTS
% AllData           output struture from segmentation 
% filter_field      a string, referring to a measured field in AllData
% filter_function   a MATLAB function (anonymous or standard) - must return a boolean
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

assert(islogical(filter_function(1)),'Filter function must be a MATLAB function that returns a logical value')


all_cond = fieldnames(AllData);
for i = 1:length(all_cond)
    if ~isfield(AllData.(all_cond{i}).Measurements,filter_field)
        error(['Error: Measurement field ''', filter_field,''' does not exist in one or more conditions.'])
    end
    filtered_cells = filter_function(AllData.(all_cond{i}).Measurements.(filter_field));
    disp([all_cond{i},': filtered out ', num2str(sum(filtered_cells)),'/', num2str(length(filtered_cells)),' cells'...,
        ' (',num2str(sum(filtered_cells)/length(filtered_cells)*100),'%)'])
    measure_names = fieldnames(AllData.(all_cond{i}).Measurements);
    AllData.(all_cond{i}).CellData(filtered_cells,:) = [];
    for j = 1:length(measure_names)
        AllData.(all_cond{i}).Measurements.(measure_names{j})(filtered_cells) = [];
    end
end

AllData_out = AllData;
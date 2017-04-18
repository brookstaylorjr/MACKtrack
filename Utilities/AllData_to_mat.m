function [fixedcell_mat, fixedcell_cols, fixedcell_cond] = AllData_to_mat(AllData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [fixedcell_mat, fixedcell_cols, fixedcell_cond] = AllData_to_mat(AllData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ALLDATA_TO_MAT converts to an (R-friendly) version of the AllData structure made by MACKtrack's
% screening functionality.
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


fixedcell_cond = fieldnames(AllData);
fixedcell_mat = [];
fixedcell_cols = {'Condition', 'Well No.','Image No.','Cell Index'};
for i = 1:length(fixedcell_cond)
    tmp_mat = AllData.(fixedcell_cond{i}).CellData(:,1:3);
    tmp_mat = [i*ones(size(tmp_mat,1),1),tmp_mat];
    
    measure_names = fieldnames(AllData.(fixedcell_cond{i}).Measurements);
    for j = 1:length(measure_names)
        var1 = AllData.(fixedcell_cond{i}).Measurements.(measure_names{j});
        if (size(var1,1) == size(tmp_mat,1) ) && (size(var1,2)==1)
            tmp_mat = [tmp_mat,var1];
            if i==1
                fixedcell_cols = cat(2,fixedcell_cols,measure_names{j});
            end
        end
    end
    fixedcell_mat = cat(1,fixedcell_mat,tmp_mat);
end

fixedcell_cond = fixedcell_cond';
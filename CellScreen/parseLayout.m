function [cond_names, cond_wells, cond_well_num] = parseLayout(layout_dir)
% [cond_names, cond_wells] = parseLayout(layout_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PARSELAYOUT looks in a specified directory for a layout.xlsx file, specifying experimental conditions in a 96-well 
% plate. It returns a list of all unique conditions, and their corresponding wells (e.g. {'A01', 'A02'}) 
% 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
% Read in spreadsheet, ignoring letter/number row & col headers
[~,~,raw] = xlsread([layout_dir,filesep,'layout.xlsx']);
raw(1,:) = []; raw(:,1) = [];
% Check remainder of spreadsheet - ensure dimensions are correct.
if (size(raw,1) ~=8) || (size(raw,2)~=12)
    error('96 well layout file is incorrect size: should have rows A-H, and columns 1-12')
end


% Blank wells will be read as 'NaN' - grab/group all conditions
cond_names = {};
cond_wells = {};
ltr = 'ABCDEFGH';
for j = 1:length(ltr)
    for k = 1:12
        if ~isnan(raw{j,k})
            tmp_val = raw{j,k};
            if isnumeric(tmp_val)
                tmp_val = num2str(tmp_val);
            end
            test_name = genvarname(tmp_val);
            % Check if condition exists (add to list if it doesn't)
            if ~ismember(test_name,cond_names)
                cond_names = cat(1,cond_names,test_name);
                cond_wells = cat(1,cond_wells,{[ltr(j),numseq(k,2)]});
            else % Condition exists - add 
                idx = strcmp(cond_names,test_name);
                cond_wells{idx} = cat(2,cond_wells{idx},{[ltr(j),numseq(k,2)]});
            end
        end
    end
end  
disp(['- - - - - - Found ''layout.xls'', under ''',layout_dir,'''. Conditions/wells: - - - - - '])
for i = 1:length(cond_names)
    if ~iscell(cond_wells{i})
        cond_wells{i} = {cond_wells{i}};
    end
    disp(['''',cond_names{i},''': from wells ' strjoin(cond_wells{i},', ')])
end    
disp('- - - - - - - - - - - - - - ')

% Optional output: specify wells as numerical array instead (i.e. 'A01' -> 101)
cond_well_num = cell(size(cond_wells));
for i = 1:length(cond_wells)
    cond_well_num{i} = [];
    for j = 1:length(cond_wells{i})
        num1 = double(cond_wells{i}{j}(1)-64);
        val1 = eval([num2str(num1),cond_wells{i}{j}(2:end)]);
        cond_well_num{i} = cat(2,cond_well_num{i},val1);
    end
end
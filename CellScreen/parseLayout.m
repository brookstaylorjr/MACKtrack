function [cond_names, cond_wells] = parseLayout(layout_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [cond_names, cond_wells] = parseLayout(layout_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PARSELAYOUT looks in a specified directory for a layout.xlsx file, specifying experimental conditions in a 96-well 
% plate. It returns a list of all unique conditions, and their corresponding wells (e.g. {'A01', 'A02'}) 
% 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
% Read in spreadsheet, ignoring letter/number row & col headers
[~,~,raw] = xlsread([layout_dir,filesep,'layout.xlsx']);
raw(1,:) = []; raw(:,1) = [];
% Blank wells will be read as 'NaN' - grab/group all conditions
cond_names = {};
cond_wells = {};
ltr = 'ABCDEFGH';
for j = 1:length(ltr)
    for k = 1:12
        if ~isnan(raw{j,k})
            test_name = genvarname(raw{j,k});
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
disp(['- - - - - - Measuring following conditions from images under directory ''',layout_dir,'''- - - - - '])
for i = 1:length(cond_names)
    disp(['''',cond_names{i},''': from wells ' strjoin(cond_wells{i},', ')])
end    
disp('- - - - - - - - - - - - - - ')
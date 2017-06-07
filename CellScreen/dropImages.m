function AllData_out = dropImages(AllData_in, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% data_out = dropImages(AllData_in, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% DROPIMAGES filters "bad" images (by string match, e.g. 'E08_s1') from an experimental set, dropping respective 
% single-cell measurements from those images. (Non-unique matches are allowed; all matched images will be displayed).
% 
% (String matches are case-sensitive)
%
% INPUTS:
% AllData_in      AllData structure (output by MACKtrack's screen analysis)
% varargin        string(s) of image(s) to filter out of AllData_in - any measurement made in a matched image will be
% dropped
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%

AllData_out = AllData_in;

% Cycle through full list of (nuclear) images - try to find (unique) match for string.
for i = 1:length(varargin)
     if ~ischar(varargin{i})
        error('Specify images using strings (e.g. ''E08_s1''.')
     end
    condition_names = fieldnames(AllData_out);
    for j = 1:length(condition_names)
        img_list = AllData_out.(condition_names{j}).Images;
        img_list = img_list(:,1);
        match_idx = find(~cellfun(@isempty,strfind(img_list,varargin{i})));
        if ~isempty(match_idx)
            s1 = ['''',varargin{i},''' found in condition ''',condition_names{j}, ''' - dropping cells from image(s):\n'];
            for k = 1:length(match_idx)
                s1 = [s1,'     - ', img_list{match_idx(k)},'\n'];
            end
            % Get indicies of cells to drop 
            drop_cells = ismember(AllData_out.(condition_names{j}).CellData(:,2),match_idx);
            % Drop these cells from measurements
            measure_names = fieldnames(AllData_out.(condition_names{j}).Measurements);
            for k = 1:length(measure_names)
                if length(AllData_out.(condition_names{j}).Measurements.(measure_names{k})) == length(drop_cells)
                    AllData_out.(condition_names{j}).Measurements.(measure_names{k})(drop_cells,:) = [];
                end
            end
            % Drop cells from CellData
            AllData_out.(condition_names{j}).CellData(drop_cells,:) = [];
            
            % Display status
            fprintf(s1)
        end    
    end
end
disp('- - - - - - - - - - - - - - -')
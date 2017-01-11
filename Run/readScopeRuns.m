function [data] = readScopeRuns(url, run_rows)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [data] = readScopeRuns(url, varargin)
%
% READSCOPERUNS looks up a published Google Sheet (by URL/ID) and converts all info in 
% the main table into a cell array, then pulls out data corresponding to specified ID(s)
%
% url         input Google Sheet URL
% run_rows    IDs corresponding to rows on Google Sheet (1 row per experimental condition)
%
% data        information from row(s) corresponding to input data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<2
    error('ERROR: specify at least 1 ID number from ''Scope Runs'' spreadsheet')
end


% Read url's contents into string
pageString = urlread(url);

% Find all tables on page.
[tables] = regexp(pageString, '(<table[^>]*>(?:(?>[^<]+)|<(?!table[^>]*>))*?</table>)','tokens');

% Build cell aray of table data

for i = 1:length(tables)
    table = tables{i};
    rows = regexpi(table, '<tr.*?>(.*?)</tr>','tokens');
    table_data = cell(0);
    % Pull off headers (if present)
    headers = regexpi(rows{1}{1}, '<th.*?>(.*?)</th>','tokens');
    if isempty(headers{1})
        start_mod = 0;
    else
        start_mod = 1;
    end
    
    % Cycle rows, then columns - pull information into a cell array
    for j = 1:(numel(rows{1})-start_mod)
        cols = regexpi(rows{1}{j+start_mod}, '<td.*?>(.*?)</td>','tokens');
        for k = 1:numel(cols{1})
            tmp = regexprep(cols{1}{k},'<.*?>', '');
            table_data{j,k} = tmp{1};
        end
    end
    
    % Stop scanning tables once we get to a valid one
    if strcmp(table_data{1,1},'#')
        break
    end
    
end
    

% Drop any empty rows/columns
table_data(:,sum(cellfun(@isempty,table_data),1)==size(table_data,1)) = [];
table_data(sum(cellfun(@isempty,table_data),2)==size(table_data,2),:) = [];


% Get IDs and pull out corresponding rows
ids = cellfun(@str2num,table_data(2:end,strcmpi(table_data(1,:),'#')));
if isempty(ids)
    error('ERROR: Couldn''t find a column named "#" in spreadsheet.')
end

[~,locs] = ismember(ids, run_rows);
locs = find(locs);

data.save_folder = table_data(2:end,strcmpi(table_data(1,:),'folder name'));
data.image_paths = table_data(2:end,strcmpi(table_data(1,:),'image path'));
data.xy_ranges = table_data(2:end,strcmpi(table_data(1,:),'xy'));
data.time_ranges = table_data(2:end,strcmpi(table_data(1,:),'t'));
data.parameter_files = table_data(2:end,strcmpi(table_data(1,:),'params file'));
data.save_dir = table_data(2:end,strcmpi(table_data(1,:),'save path'));
data.modify = table_data(2:end,strcmpi(table_data(1,:),'other params'));

% Double check that all necessary data was found and subset data
types = fieldnames(data);
for i = 1:length(types)
    if isempty(data.(types{i}))
        error(['ERROR: Couldn''t find column "',types{i},'" in spreadsheet.'])
    end
    data.(types{i}) = data.(types{i})(locs);
end

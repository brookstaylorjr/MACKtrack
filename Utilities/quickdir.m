function file_list = quickdir(directory)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% file_list = quickdir(directory)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% QUICKDIR uses a system directory command to quickly get file contents of a target directory. Filters hidden files
% and system shortcuts (e.g. 'Thumbs.db', '.' and '..')
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%

if ~exist(directory,'dir')
    error('QUICKDIR requires a valid directory')
end


% Windows-specific
if ispc
    [~, filenames] = system(['dir /b "',directory,'"']);
else % OSX/Linux
    [~, filenames] = system(['ls -1 "', directory,'"']);
end


file_list = strsplit(filenames,'\n')';
filter1 = @(a) isempty(a) || strcmp(a(1),'.') || strcmp(a,'Thumbs.db')|| strcmp(a,'desktop.ini');
file_list(cellfun(filter1,file_list)) = [];

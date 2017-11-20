function file_list = quickdir(directory,quick_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% file_list = quickdir(directory)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% QUICKDIR uses a system directory command to quickly get file contents of a target directory. Filters hidden files
% and system shortcuts (e.g. 'Thumbs.db', '.' and '..')
%
% If quick_flag is enabled (off by default), command will timeout to an empty string (OSX/Linux only). (Uses Perl)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%

if nargin<2
    quick_flag = 0;
end

if ~exist(directory,'dir')
    error('QUICKDIR requires a valid directory')
end


% Windows-specific
if ispc
    [~, filenames] = system(['dir /b "',directory,'" /O:N']);
else % OSX/Linux
    if ~quick_flag
        [~, filenames] = system(['ls -1 "', directory,'"']);
    else
        [status, filenames] = system(['perl -e "alarm 20; exec @ARGV" "ls -1 "',directory,'""']);
        if status ~=0
            filenames = 'note: too many files to show! Use ''...'' to navigate further';
        end
        
    end
end


file_list = strsplit(filenames,'\n')';
filter1 = @(a) isempty(a) || strcmp(a(1),'.') || strcmp(a,'Thumbs.db')|| strcmp(a,'desktop.ini');
file_list(cellfun(filter1,file_list)) = [];

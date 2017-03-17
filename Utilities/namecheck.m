function name_out = namecheck(name_in,default_empty)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% name_out = namecheck(name_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NAMECHECK   provides system-specific formatting of a file name/path. If input is blank,
% will return a fallback path (default fallback is current working directory).
%
% INPUT:
% name_in         input filename
% default_empty   output if an empty string is passed (default is current directory)
%
% OUTPUT:
% name_out   output (cleaned) filename
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<2
    default_empty = [pwd,filesep];
end

if isempty(name_in)
    name_out = default_empty;
    return;
end

% Convert filesep characters
os = computer;
if strcmp(os(1:4),'MACI') || strcmp(os(1:4),'GLNX') % *NIX - convert \ to /
    if ~isempty(strfind(name_in,'\'))
        disp('Note: converting ''\'' to ''/'' in file path')
        name_in(strfind(name_in,'\')) = filesep;
    end
    
else % Windows - convert / to \
    if ~isempty(strfind(name_in,'/'))
        name_in(strfind(name_in,'/')) = filesep;
    end
end

% Convert all double slashes in filename to single slashes (exclude 1st, in case it was a network path of some kind)
idx = strfind(name_in,[filesep,filesep]);
idx(idx==1) = [];
while ~isempty(idx)
    name_in(idx) = '';
    idx = strfind(name_in,[filesep,filesep]);
    idx(idx==1) = [];
end


% Add filesep to end of a directory, convert a filename to an absolute path
if exist(name_in,'dir')
    if ~strcmp(name_in(end),filesep)
        name_in = [name_in,filesep];
    end
elseif isempty(strfind(name_in,filesep))
    name_in = which(name_in);
end


name_out = name_in;
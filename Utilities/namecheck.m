function name_out = namecheck(name_in,default_empty)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% name_out = namecheck(name_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NAMECHECK   provides system-specific checking of a file name/path.
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
        warn('Note: converting ''\'' to ''/'' in file path')
        name_in(strfind(name_in,'\')) = filesep;
    end
    
else % Windows - convert / to \
    if ~isempty(strfind(name_in,'/'))
        name_in(strfind(name_in,'/')) = filesep;
    end
    
end

% Check if file exists; convert to absolute path
if ~exist(name_in,'file') && ~exist(name_in,'dir')
    error(['Couldn''t find file/directory ''',name_in,''' - exiting.'])
end


% Add filesep to end of a directory
if exist(name_in,'dir')
    if ~strcmp(name_in(end),filesep)
        name_in = [name_in,filesep];
    end
elseif isempty(strfind(name_in,filesep))
    name_in = which(name_in);
end


name_out = name_in;
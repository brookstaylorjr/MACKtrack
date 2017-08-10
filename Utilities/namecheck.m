function name_out = namecheck(name_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% name_out = namecheck(name_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NAMECHECK   provides system-specific formatting of a file name/path.
%
% INPUT:
% name_in         input filename
% default_empty   output if an empty string is passed (default is current directory)
%
% OUTPUT:
% name_out   output (cleaned) filename
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Convert filesep characters
os = computer;
if strcmp(os(1:4),'MACI') || strcmp(os(1:4),'GLNX') % *NIX - convert \ to /
    if ~isempty(strfind(name_in,'\'))
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


% Add filesep to end of a directory
if exist(name_in,'dir')
    if ~strcmp(name_in(end),filesep)
        name_in = [name_in,filesep];
    end
end


name_out = name_in;
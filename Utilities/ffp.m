function name = ffp(name)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FFP (format file path) converts a file path slashes to either windows or
% Linux/Unix. Duplicate slashes are also formatted out.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if ispc
    name(name=='/') = '\';
else
    name(name=='\') = '/';
end
idx = strfind(name,[filesep filesep]);
while ~isempty(idx)
    name(idx) = '';
    idx = strfind(name,[filesep filesep]);
end
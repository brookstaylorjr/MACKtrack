function combined_struct = combinestructures(varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% combined_struct = combinestructures(struct1, struct2,....)
%
% COMBINESTRUCTURES combines the fields of multiple structures (they MUST be 1x1). Each 
% structure is passed as a separate argument. If there is field overlap, function will
% generate a warning, and use the earlier-listed structure's field.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

for i = 1:length(varargin)
    if length(varargin{i}) > 1
       error('Error in ''combinestructures'': structures may not have a size > 1')
    end
end

combined_struct = varargin{1};

if nargin>1
    for i = 2:length(varargin)
        names = fieldnames(varargin{i});
        for j = 1:length(names)
           if ~isfield(combined_struct,names{j})
               combined_struct.(names{j}) = varargin{i}.(names{j});
           elseif ~isequal(combined_struct.(names{j}),varargin{i}.(names{j}))
               warning(['Warning (''combinestructures''): structures have same field (''',names{j},''') / different data. Using 1st struct''s info'])         
           end
        end
    end
end
function [names_out, scope_type] = wellmatch(contents, well_in, channel_in, scope_type)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% names_out = wellmatch(names_in, well_in, channel_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% WELLMATCH performs partial name matching on a (cell matrix) list of file contents, returning the files that match
% the provided well and channel.
%
% Valid naming schema:
% 
% (1): MetaXpress (Molecular Devices) ('metaexpress')
%   A) .tif only
%   B) anything with the word "thumb" in the title will be ignored
%   C) wells are zero-padded and surrounded by underscores (e.g. "_H08_")
%   D) filenames include a database ID - a set of 8,4,4,4, and 12 hex characters, separated by dashes
%
% (2): Slidebook (3i) ('slidebook')
%   A) .tif or .tiff
%   B) wells are NOT zero-padded, are preceded by a space, and succeeed by an underscore (e.g. " H5_" 
%   C) omit '.log' and '.xml' files (if there)
%
% (3): Steve Cappell's naming script ('cappell')
%   A) .tif only
%   B) wells are specified by (non-zero-padded) numbers @ the beginning - e.g. 5_3_1_DAPI.tif is E03, site 5.
%
% (4): Nikon NIS elements HCA pipeline ('nikon')
%   A) tif only
%   B) wells are zero-padded (e.g. 'H04')
%   C) contains term (e.g.) 'WellA01'
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%%
% 1) Determine what kind of images we're dealing with (if not provided)
if nargin < 4
    imglist  = contents(~cellfun(@isempty,strfind(contents,'.tif')));
    imglist  = imglist(cellfun(@isempty,strfind(imglist,'.xml')));
    imglist  = imglist(cellfun(@isempty,strfind(imglist,'.log')));
    imglist = imglist(cellfun(@isempty,strfind(imglist,'_thumb')));

    % a) MetaXpress: image names and should have database ID @ end - a set of 8,4,4,4, and 12 characters
    hyphen_places = diff([strfind(imglist{1},'-'),strfind(imglist{1},'.tif')]);
    if (length(hyphen_places)>=4) && isequal(hyphen_places(end-3:end),[5 5 5 13])
        scope_type = 'metaxpress';
    % b) Slidebook: last characters of image names should be '_C[0-9].tif'
    elseif ~isempty(regexp(imglist{1}, '_C[0-9].tif','ONCE'))
        scope_type = 'slidebook';
    % c) Steve Cappell's script: image names should begin with single number-underscore-number
    elseif ~isempty(regexp(imglist{1}(1:3), '[1-8]_[1-9]','ONCE'))
        scope_type = 'cappell';
    % d) NIS elements: image name should begin with 'Well'
    elseif  ~isempty(regexp(imglist{1}, 'Well[A-H][0-1]','ONCE'))
        scope_type = 'nikon';
    else
        error('No valid image names found in your specified directory')
    end
end


%% 2) Based on scope_type, filter images appropriately
if ~isempty(well_in)
    switch scope_type
        case 'metaxpress'
            names_out = contents(cellfun(@isempty,strfind(contents,'_thumb'))...
                &~cellfun(@isempty,strfind(contents,['_',well_in,'_'])));     
        case 'slidebook'
            % Strip zero padding from well
            if strcmp(well_in(2),'0')
                well_in = well_in([1 3]);
            end
            names_out = contents(cellfun(@isempty,strfind(contents,'.xml'))...
                & cellfun(@isempty,strfind(contents,'.log'))...
                &~cellfun(@isempty,strfind(contents,[' ',well_in,'_'])));
        case 'cappell'
            % Convert well_in to numerical value
            row = double(well_in(1))-64;
            col = eval(well_in(2:end));
            well_in = [num2str(row),'_',num2str(col)];
            firstfour = @(str) str(1:4);
            contents_tmp = cellfun(firstfour,contents,'UniformOutput',0);        
            names_out = contents(~cellfun(@isempty,strfind(contents_tmp,well_in)));
        case 'nikon'
            names_out = contents(~cellfun(@isempty,strfind(contents,['Well',well_in,'_'])));
    end
else % No well specified - just filter out 'bad' files
    switch scope_type
        case 'metaxpress'
            names_out = contents(cellfun(@isempty,strfind(contents,'_thumb')));     
        case 'slidebook'
            names_out = contents(cellfun(@isempty,strfind(contents,'.xml'))...
                & cellfun(@isempty,strfind(contents,'.log')));
        otherwise     
            names_out = contents;   
    end
    

end


% If input channel is defined, filter further
if nargin>2
    names_out = names_out(~cellfun(@isempty,strfind(names_out,channel_in)));
end

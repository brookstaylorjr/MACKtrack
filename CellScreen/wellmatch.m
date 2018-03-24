function [names_out, scope_type] = wellmatch(contents, well_in, channel_in, scope_type)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% names_out = wellmatch(names_in, well_in, channel_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% WELLMATCH performs partial name matching on a (cell matrix) list of file contents, returning the files that match
% the provided well and channel.
%
% Valid schema, current as of 9/13/2017:
% 
% Schema 1: MetaXpress (Molecular Devices)
% 1) .tif only
% 2) anything with the word "thumb" in the title will be ignored
% 3) wells are zero-padded and surrounded by underscores (e.g. "_H08_")
%
%
% Schema 2: Slidebook (3i)
% 1) .tif only (NOT compatible w/ old versions of Slidebook!!)
% 2) wells are NOT zero-padded, are preceded by a space, and succeeed by an underscore (e.g. " H5_" 
% 3) omit '.log' and '.xml' files (if there)
%
% Schema3: Steve Cappell's naming script
% 1) .tif only
% 2) wells are specified by (non-zero-padded) numbers - e.g. 5_03_01_DAPI.tif is E03, site 5.
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% 1) Determine what kind of images we're dealing with
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
    else
        error('No valid image names found in your specified directory')
    end


end


%% 2) Based on scope_type, filter images appropriately
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
        
end

% If input channel is defined, filter further
if nargin>2
    names_out = names_out(~cellfun(@isempty,strfind(names_out,channel_in)));
end

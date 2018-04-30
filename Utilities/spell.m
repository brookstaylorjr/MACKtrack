function [str] = spell(n,precision)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [str] = spell(n,precision)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SPELL formats and displays a matrix in a clean, copyable format. SPELL uses precision of 3 by default; can pass as an 
% additional (optional) argument (default precision = 4)
%
% INPUTS:
% n           2D matrix of numerical values
% precision   (optional) numerical precision of output values (default=4)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin==1
    precision = 3;
end

if isempty(n)
    str = '';
   
else
    tot_length = 3+precision;

    if size(n,1)>1
        disp('[...')
        for i=1:size(n,1)
            str = [];
            for j = 1:size(n,2)
                substr = num2str(n(i,j),precision);
                str = [str,substr,repmat(' ',1,tot_length-length(substr))];
            end
            disp(str)
        end
        disp('];')
    else
        str = [];
        for j = 1:size(n,2)
            substr = num2str(n(1,j),precision);
            str = [str,substr,repmat(' ',1,tot_length-length(substr))];
        end


    end
end

        if nargout<1
            disp(['[ ', str,' ]']);
        else
            str = ['[ ', str,' ]'];
        end
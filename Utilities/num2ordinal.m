function chr = num2ordinal(num,isd)
% Convert a numeric array to a character array of integers and ordinal suffixes (string).
%
% (c) 2016 Stephen Cobeldick
%
% Convert a numeric array into a character array containing integers and
% ordinal suffixes, e.g. 1 -> '1st'. Fully vectorized for operating on arrays!
%
% Syntax:
% chr = num2ordinal(nun)
% chr = num2ordinal(num,isd)
%
% NUM2ORDINAL rounds any fractional values to the nearest integers.
% It also provides the correct suffixes for values ending in 11, 12 or 13.
%
% See also NUM2WORDS WORDS2NUM SPRINTF INT2STR NUM2STR NUM2CELL NATSORT NUM2BIP NUM2SIP CELLSTR STRTRIM STRCAT INTMAX
%
% ### Examples ###
%
% num2ordinal(1)
%  ans = '1st'
%
% num2ordinal(1:6)
%  ans = ['1st';'2nd';'3rd';'4th';'5th';'6th']
%
% num2ordinal([1,11,111,1111])
%  ans = ['   1st';'  11th';' 111th';'1111th']
%
% num2ordinal(100:113,false)
%  ans = ['th';'st';'nd';'rd';'th';'th';'th';'th';'th';'th';'th';'th';'th';'th']
%
% num2ordinal(intmax('int64')-4)
%  ans = '9223372036854775803rd'
%
% num2ordinal([-0,0])
%  ans = ['-0th';' 0th']
%
% vals = [-1,-0,0;-Inf,NaN,Inf];
% reshape(strtrim(cellstr(num2ordinal(vals))),size(vals))
%  ans = {'-1th','-0th','0th';'-Infth','NaNth','Infth'}
%
% ### Input & Output Arguments ###
%
% Inputs (*==default):
% num = NumericArray (any size), with values to convert to string.
% isd = LogicalScalar, *true/false -> digits+suffix / suffix only.
%
% Output:
% chr = CharArray, integers & ordinal suffixes of the rounded <num> values,
%       where the rows of <chr> are linear indexed from <num>.
%
% chr = num2ordinal(num,*isd)

typ = class(num);
switch typ(1:3)
    case {'dou','sin'}
        assert(isreal(num),'First input <num> must be a real numeric.')
        fmt = '.0f';
    case 'uin'
        fmt = 'lu';
    case 'int'
        fmt = 'ld';
    otherwise
        error('First input <num> must be a numeric scalar/vector/matrix/array.')
end
num = num(:);
%
% Convert from number to character array plus 'th' suffix:
if isempty(num)
    % Empty numeric.
    chr = '';
    return
    %
elseif nargin>1&&islogical(isd)&&~isd
    % Return suffixes only.
    len = 2;
    chr(1:numel(num),2) = 'h';
    chr(:,1) = 't';
    %
elseif isscalar(num)
    % Numeric scalar.
    chr = sprintf(['%',fmt,'th'], num);
    len = numel(chr);
    %
else
    % Numeric matrix.
    % Calculate the number of digits required:
    len = 1 + max(0, floor(log10(double(abs(num)))));
    len(~isfinite(num)) = 3;
    len = 2 + max(len + (num<0 | 1./num<0));
    % Convert from numeric to character:
    fmt = sprintf('%%%.0f%sth', len-2, fmt);
    chr = reshape(sprintf(fmt, num), len, []).';
    %
end
%
% Determine the correct ordinal suffixes:
rm = rem(num,100);
th = 10.5<=rm & rm<13.5;
rm = rem(num,10);
st = ~th & 0.5<=rm & rm<1.5;
nd = ~th & 1.5<=rm & rm<2.5;
rd = ~th & 2.5<=rm & rm<3.5;
%
% Replace 'th' suffixes with the correct suffixes:
leu = len-1;
chr(st,leu) = 's';
chr(st,len) = 't';
chr(nd,leu) = 'n';
chr(nd,len) = 'd';
chr(rd,leu) = 'r';
chr(rd,len) = 'd';
%
end
%----------------------------------------------------------------------
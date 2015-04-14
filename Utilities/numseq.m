function [numstring] = numseq(num1, digits)
% NUMSEQ Converts numerical value num1 (e.g. 1) to a string of length 'digits'
% Example: numstring = numseq(15, 4) returns
%       numstring = 
%       '0015'
%
% Inputs
% -num1: input number
% -digits: total length of output string (pads zeros on front)
%
% Outputs
% -numstring: output string made from num1


numstring = '';
lengthNum1 = length(num2str(num1));
for i = 1:(digits-lengthNum1)
numstring= strcat(numstring,'0');
end
numstring = strcat(numstring,num2str(num1));

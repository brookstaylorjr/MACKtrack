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

base_len = length(num2str(num1));
numstring = [repmat('0',[1,digits-base_len]),num2str(num1)];
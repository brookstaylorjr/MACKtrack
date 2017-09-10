function [numstring] = numseq(num_in, len)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [numstring] = numseq(num_in, len)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NUMSEQ Converts numerical value num_in to a string of length 'digits'
% Note: for decimal values, values after decimal will NOT be considered in overall string
% length.
%
% Example 1 : numstring = numseq(15, 4) returns '0015'
%
% Example 1 : numstring = numseq(15.01, 4) returns '0015.01'
%
%
% INPUTS
% num_in     input number
% len        total length of output string (pads zeros on front)
%
% OUTPUTS
% numstring  output string made from num1
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

base_len = length(num2str(floor(num_in)));
numstring = [repmat('0',[1,len-base_len]),num2str(num_in)];
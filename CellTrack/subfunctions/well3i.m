function string_out = well3i(num_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% string_out = well3i(num_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% WELL3i name converts between a row/col/site value and a valid image name (as 
% created/exported by Slidebook's multiwell + export funcionality)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% Define row (A-H)
row = floor(num_in/100);

% Define column (1-12, zero-padded)
col = rem(floor(num_in),100);

% Define site (1-XX, zero-padded to allow up to 99 sites)
site = round(rem(num_in,1)*100);


% Error checking on num_in. Criteria:
% 1) Hundreds digit is btw 1-8
% 2) Tens/ones digits are btw 1-12
% 3) Sites are btw 1-99
is_valid = @(row,col,site) (row>=1)&&(row<=8)&&(col>=1)&&(col<=12)&&(site>=1)&&(site<=99);

if ~is_valid(row,col,site)
    disp('- - - - - - - - - - ')
    disp(['Note: site [',num2str(num_in),'] is specified incorrectly. Should be specified as: [XYY.ZZ], where'])
    disp('    - X is well row (converted to number 1-8)')
    disp('    - YY is well column')
    disp('    - ZZ is the site number]')
    disp('    e.g.: Well A8, site 5 -> 108.05')
    disp('- - - - - - - - - - ')

end


% Assemble string out (subdirectory + image prefix)
string_out = [char(64+row),num2str(col),'_',numseq(site-1,3)]; % 3i sites start @ 0, not 1

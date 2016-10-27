function [] = printfigs(savedir, savestem, figs, figs2)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = printfigs(savedir, savestem)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PRINTFIGS scans two structures of figure hangle objects (figs and figs2) - each field refers to a unique and open
% figure. All figures in 'figs' will be printed as .eps files, and all figures in 'figs2' will be printed as a .png.
%
% INPUTS
% savedir     directory to save figures to
% savestem    naming stem for each figure
% figs 
% figs2
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Checks:
if nargin<4
    figs2 = struct;
end

%% Ensure save directory exists, then append filesep to end (if not already present)
if ~exist(savedir,'dir')
    error(['Specified save directory, ',savedir, 'doesn''t exist'])
else
    if ~strcmp(savedir(end),filesep)
        savedir = [savedir,filesep'];
    end
    disp(['Saving all figs to ', savedir,'...'])
end

% Save figs without transparency as .eps
subfigs = fieldnames(figs);
for i = 1:length(subfigs)
    for j = 1:numel(figs.(subfigs{i}))  
        figname = [savestem,'_',subfigs{i},'(',num2str(j),').eps'];
        print(figs.(subfigs{i})(j),[savedir,figname], '-depsc')      
        disp([' - printed "',figname,'"'])
    end
end
% Save figs with transparency as .png (@ high resolution)
subfigs = fieldnames(figs2);
for i = 1:length(subfigs)
    for j = 1:numel(figs2.(subfigs{i}))
        figname = [savestem,'_',subfigs{i},'(',num2str(j),').svg'];
        print(figs2.(subfigs{i})(j),[savedir,figname], '-dsvg','-painters')
        disp([' - printed "',figname,'"'])
    end
end
disp('- - - - - - - - - - - - - - - - - - - ')
%%
global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

dirs = frGetDirs;
addpath(fullfile(dirs.svn,'\third\CellSort 1.1'));

%% multi-plane registration registration
newdir = '2010\110304\LG26_OriA';
newdir = '2010\100923\OriA';

expt = frGetExpt(newdir);
frRegister2(newdir,'Overwrite',true,'Oversampling',10,...
               'DoRecurse',false,'Engine','subpixel','nPlanes',8);


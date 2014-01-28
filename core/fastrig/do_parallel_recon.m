function do_parallel_recon(exptName, seriesName)

rootPath = 'i:/users/histed/matlab/reidlab-full-svn/users/histed/py-work';
ncpus = 4;

% add root dir to path to find DLLs

%tStr = sprintf('!%s', fullfile(rootPath, 'do_recon_single.py'));
%eval(tStr);

commandStr = sprintf('%s %s %s %d &', ...
                     fullfile(rootPath, 'do_recon_single.py'), ...
                     exptName, ...
                     seriesName, ...
                     ncpus);

dos(commandStr);

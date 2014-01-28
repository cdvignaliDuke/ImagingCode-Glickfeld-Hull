%% consts
exptName = 'mouse090619';
frameTimeMs = frGetFrameRate;
seriesName = 'vstim2';

doOverwrite = 0;
dirs = frGetDirs;  % you may need to change dirs; see frGetDirs.m
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ... 
    'DoRecurse', false, ...
    'WaitForDads', true, ...  % if scan is aborted, set false
    'AverageN', 1, ...
    'RecalculateDelay', false);

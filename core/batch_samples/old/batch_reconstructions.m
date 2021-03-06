
%% consts
exptName = 'mouse090619';
frameTimeMs = frGetFrameRate;



%% do reconstruction
seriesName = 'vstim6';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
%dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim7';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim8';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim9';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim11';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim12';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim13';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim14';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim15';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim16';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim17';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim18';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

seriesName = 'vstim19';
doOverwrite = 0;
dirs = frGetDirs;
dadDir = sprintf('%s/%s', exptName, seriesName);
frReconstruct(dadDir, ...
    'Overwrite', doOverwrite, ...
    'DoRecurse', false, ...
    'WaitForDads', true, ...
    'AverageN', 1, ...
    'RecalculateDelay', false);

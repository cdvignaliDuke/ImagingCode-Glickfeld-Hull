function exptNames = frListExpts(animal)
%frListExpts (fastrig): get all data directories for a given animal
%   exptNames = frListExpts(animal)

dirs = frGetDirs;
d = dir(fullfile(dirs.data, animal));

dNames = {d.name};
notDotIx = cellfun(@isempty, regexpi('^\.\.?$', dNames, 'match'));

exptNames = dNames(notDotIx);

function dadName = makeDadName(tPath, tNum)
%$Id$
[pP pD pE] = fileparts(tPath);
dadName = fullfile(tPath, sprintf('%s_%06d.dad', pD, tNum));
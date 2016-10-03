function runOriCellSets(datasetStr,cellsOnly)
eval(['awFSAVdatasets' datasetStr])
rc = behavConstsAV;
for iexp = 1:length(expt)
    [cellsSelect, OSI, DSI] = OriCellSets(rc, expt, iexp,cellsOnly);
end
end
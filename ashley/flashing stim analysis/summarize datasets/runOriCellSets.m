function runOriCellSets(datasetStr,cellsOnly)
if isempty(datasetStr)
    eval(['awFSAVdatasets' datasetStr])
    isFSAV = 0;
elseif strcmp(datasetStr(1:3),'FSA')
    eval(datasetStr)
    isFSAV = 1;
else
    eval(['awFSAVdatasets' datasetStr])
    isFSAV = 0;
end
rc = behavConstsAV;
for iexp = 1:length(expt)
    [cellsSelect, OSI, DSI] = OriCellSets(rc, expt, iexp,cellsOnly,isFSAV);
end
end
function getOriTuningFSAV(ds)
rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
eval(dataGroup)

tuningDownSampFactor = 10;
nexp = size(expt,2);

for iexp = 1:nexp
    disp(iexp)
    dataPath = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,...
        'two-photon imaging', expt(iexp).date, expt(iexp).dirtuning);
    load(fullfile(dataPath,'timecourses.mat'))
    mworks = loadMworksFile(...
        expt(iexp).SubNum,expt(iexp).date,expt(iexp).dirtuning_time);
    
    [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability] = ...
    getOriTuning(data_tc_subnp,mworks,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));
    save(fullfile(dataPath,'oriTuningAndFits'),...
        'avgResponseEaOri','semResponseEaOri',...
        'vonMisesFitAllCells','fitReliability')
end
end
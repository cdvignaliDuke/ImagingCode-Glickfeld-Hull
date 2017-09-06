awFSAVdatasets_V1
rc = behavConstsAV;
master_file_AW0

for iexp = 1:size(expt,2)
    subNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    disp([mouse ' ' expDate])
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];

        % load cell timecourses and mworks file
        fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
        load(fullfile(fn,'timecourses.mat'))

        % combine data
        data = dataTimecourse.dataTCsub;
        clear dataTimecourse
        [nfr,nc] = size(data);

        % get spike times
        spikeTimes = deconvolution_standalone2(ops0,data);
        save(fullfile(fn,'spikeTimes'),'spikeTimes')
    end
end
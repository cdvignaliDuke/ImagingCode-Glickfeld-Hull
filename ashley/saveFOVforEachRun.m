rc = behavConstsAV;
awData_audMod_V13trialtypes
slct_expt = 3;

nframes = 1000;
for iexp = slct_expt
    subnum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    fn = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        fnrun = fullfile(fn,runFolder);

        [~,data] = Load_SBXdataPlusMWorksData(...
            subnum,expDate,runTime,mouse,runFolder,fName,nframes);

        load(fullfile(fn,'regOuts&Img'));
        [~,dataReg] = stackRegister(data,data_corr_img);

        FOV = mean(dataReg,3);
        writetiff(FOV,fullfile(fnrun,'FOV'));
    end
end


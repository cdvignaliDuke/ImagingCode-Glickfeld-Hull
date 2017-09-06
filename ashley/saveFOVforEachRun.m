rc = behavConstsAV;
awFSAVdatasets_V1gad
slct_expt = 3;

nframes = 1000;
for iexp = slct_expt
    subnum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    fn = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);
    
    offset = 0;
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        runTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        fnrun = fullfile(fn,runFolder);

        [data,nframestotal] = loadsbx_choosepmt(1,...
            mouse,expDate,runFolder,fName,nframes);
        
        if expt(iexp).greenredsimultaneous
            load(fullfile(fn,'reg2RedOuts&Img'));
            outs = double(out_bx(1+offset:nframestotal+offset,:));
            [~,dataReg] = stackRegister_MA(data,[],[],outs);
        else
            load(fullfile(fn,'regOuts&Img'));
            [~,dataReg] = stackRegister(data,data_corr_img);
        end
        offset = offset+nframestotal;
        FOV = mean(dataReg,3);
        writetiff(FOV,fullfile(fnrun,'FOV'));
    end
end


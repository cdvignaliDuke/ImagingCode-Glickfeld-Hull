clear all
close all
ds = 'movDotsSpeedTun_V1';
slct_exp = [1];
%%
rc = behavConstsAV;
eval(ds)
for iexp = slct_exp

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
down = expt(iexp).downSampleRate;

for irun = 1:expt(iexp).nrun
    runFolder = expt(iexp).runs(irun,:);
    runTime = expt(iexp).time_mat(irun,:);
    fnrun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,runFolder);
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);
    %% load tuning data
    fName = [runFolder '_000_000'];
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,runTime,mouse,runFolder,fName);  

    % down-sample
    data_down = stackGroupProject(data,down);
    clear data

    % remove negative data by subtraction
    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    clear data_down

    %% load outs and re-regester data
    load(fullfile(fnrun,'regOuts&Img.mat'));

    [~, data_reg] = stackRegister_MA(data_sub,[],[],out_reg);

    %% load mask

    load(fullfile(fn,'final_mask.mat'))

    %% get timecourses and subtract neuropil
    buf = 4;
    np = 6;
    data_tc = stackGetTimeCourses(data_reg,mask_cell);
    data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);


    %% save data
    save(fullfile(fnrun,'timecourses.mat'),'data_tc_subnp')
    save(fullfile(fnrun,'raw_tc.mat'),'out_reg','data_tc','buf','np')
end
end
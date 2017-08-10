clear all
close all
ds = 'awData_audMod_V13trialtypes';
doRedChannelOnly = 0;
slct_exp = 1;
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
    if down > 1
        data_down = stackGroupProject(data,down);
        clear data
    else
        data_down = data;
        clear data
    end

    % remove negative data by subtraction
    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    clear data_down

    %% load outs and re-regester data
    load(fullfile(fnrun,'regOuts&Img.mat'));

    [~, data_reg] = stackRegister_MA(data_sub,[],[],out_reg);

    %% load mask
    if doRedChannelOnly
        load(fullfile(fn,'red_mask.mat'))
        mask_cell = red_mask;
    else
        try
            load(fullfile(fn,'final_mask.mat'))
        catch
            load(fullfile(fn,runFolder,'final_mask.mat'))
        end
    end
    %% get timecourses and subtract neuropil
    buf = 4;
    np = 6;
    data_tc = stackGetTimeCourses(data_reg,mask_cell);
    data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);


    %% save data
    if doRedChannelOnly
        save(fullfile(fnrun,'red_timecourses.mat'),'data_tc_subnp')
        save(fullfile(fnrun,'raw_tc.mat'),'out_reg','data_tc','buf','np')
    else
        save(fullfile(fnrun,'timecourses.mat'),'data_tc_subnp')
        save(fullfile(fnrun,'raw_tc.mat'),'out_reg','data_tc','buf','np')
    end
end
end
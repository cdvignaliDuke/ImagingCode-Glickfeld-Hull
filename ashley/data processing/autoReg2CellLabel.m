clear all
close all
ds = '_V1som';
rc = behavConstsAV;
eval(['awFSAVdatasets' ds])
slct_expt = 1;
%%
for iexp = slct_expt

    SubNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    dirFolder = expt(iexp).dirtuning;
    dirTime = expt(iexp).dirtuning_time;
    redFolder = expt(iexp).redChannelRun;


    %% load and register red channel
    redData = loadsbx_choosepmt(2,mouse,expDate,redFolder,[redFolder '_000_000']);
    redDown = stackGroupProject(redData,10);
    [~,redReg] = stackRegister(redDown,redDown(:,:,1));
    redImage = mean(redReg,3);
    %% load and register all behavior data
    % concatenate data files
    for irun = 1:expt(iexp).nrun

        runFolder = expt(iexp).runs(irun,:);
        expTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        [input, data_temp] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName);

        if irun == 1
            data_bx = data_temp;
            input_bx = input;
        else
            data_bx = cat(3,data_bx,data_temp);
            try
                input_bx = [input_bx input];
            catch
                inpNames1 = fieldnames(input_bx);
                inpNames2 = fieldnames(input);
                inpLong = gt(length(inpNames1),length(inpNames2));
                if inpLong == 1
                    inpPlusInd = ismember(inpNames1,inpNames2);
                    inpPlus = inpNames1(~inpPlusInd);
                    for i = 1:length(inpPlus)
                        input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
                    end
                else
                    inpPlusInd = ismember(inpNames2,inpNames1);
                    inpPlus = inpNames2(~inpPlusInd);
                    for i = 1:length(inpPlus)
                        input_bx.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames1));
                    end
                end
                input_temp = [input_bx input];
            end
        end
        clear data_temp input
    end
    input_bx = concatenateDataBlocks(input_bx);
    nfr_bx = size(data_bx,3);

    % remove negative data by subtraction
    data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
    data_bx = data_sub;
    clear data_sub

    % register all behavior data
    [out_bx data_bx_reg] = stackRegister(data_bx,redImage);
    clear data_bx


    % load tuning data
    down = 10;
    fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
    fName = [dirFolder '_000_000'];
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);

    % down-sample
    data_tun_down = stackGroupProject(data_tun,down);
    clear data_tun

    % remove negative data by subtraction
    data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
    data_tun_down = data_sub;
    clear data_sub

    [out_tun data_tun_reg] = stackRegister(data_tun_down,redImage);
    clear data_tun

    fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fnout,'dir')
        mkdir(fnout)
    end
    save(fullfile(fnout,'regOuts&Img.mat'),'out_bx','out_tun', 'redImage')

    %% ********* behavior max dF/F ************

    fnbx = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    % get behavior max dF/F
    taskMaxDFF

    %% ********* tuning max dF/F ************
    %%
    nfr_tun = size(data_tun_reg,3);
    fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);

    dirTuningMaxDFF
end
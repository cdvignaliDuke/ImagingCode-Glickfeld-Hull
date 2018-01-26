clear all
close all
ds = 'awFSAVdatasets_V1gad';
rc = behavConstsAV;
eval(ds)
slct_expt = 10:size(expt,2);
%%
for iexp = slct_expt
    
    SubNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    dirFolder = expt(iexp).dirtuning;
    dirTime = expt(iexp).dirtuning_time;
    redChannelOn = expt(iexp).greenredsimultaneous;
    regImgInd = expt(iexp).regImgStartFrame;
    down = 10;

    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end
    fnout = fullfile(fn,'data processing');
    %% load and register all behavior and tuning data
    % concatenate data files
    for irun = 1:expt(iexp).nrun
        runFolder = expt(iexp).runs(irun,:);
        expTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];

        if redChannelOn == 1
            data_temp = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName);
        else
            data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
        end
        if irun == 1
            data_bx = data_temp;
        else
            data_bx = cat(3,data_bx,data_temp);
        end
        clear data_temp 
    end

    nfr_bx = size(data_bx,3);
    if isnan(regImgInd)
        load(fullfile(fn,'regOuts&Img.mat'))
        regImg = data_corr_img;
    else
        regImg = mean(data_bx(:,:,regImgInd:(regImgInd+100)),3);
    end

    data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
    data_bx = data_sub;
    clear data_sub

    [out_bx,data_bx_reg] = stackRegister(data_bx,regImg);
    clear data_bx

    data_bx_down = stackGroupProject(data_bx_reg,down);
    clear data_bx_reg

    if redChannelOn == 1
        for irun = 1:expt(iexp).nrun
            runFolder = expt(iexp).runs(irun,:);
            expTime = expt(iexp).time_mat(irun,:);
            fName = [runFolder '_000_000'];

            data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
            
            if irun == 1
                data_bx = data_temp;
            else
                data_bx = cat(3,data_bx,data_temp);
            end
            clear data_temp
        end
        data_bx_down = stackGroupProject(data_bx,down);
        clear data_bx
        [~,data_bx_reg] = stackRegister(data_bx_down,regImg);
        corrImage_bx = getPixelCorrelationImage(data_bx_reg);
        clear data_bx_reg data_bx_down
        
    else
        corrImage_bx = getPixelCorrelationImage(data_bx_down);
        clear data_bx_down
    end

    % tuning data
    fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
    fName = [dirFolder '_000_000'];

    if redChannelOn == 1
        data_tun = loadsbx_choosepmt(2,mouse,expDate,dirFolder,fName);
    else
        data_tun = loadsbx_choosepmt(1,mouse,expDate,dirFolder,fName);
    end

    data_tun_down = stackGroupProject(data_tun,down);
    clear data_tun

    data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
    clear data_tun_down

    [out_tun,data_tun_reg] = stackRegister(data_sub,regImg);
    clear data_sub
    if redChannelOn == 1
        data_tun = loadsbx_choosepmt(1,mouse,expDate,dirFolder,fName);
        data_tun_down = stackGroupProject(data_tun,down);
        clear data_tun
        data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
        clear data_tun_down
        [out_tun,data_tun_reg] = stackRegister(data_sub,regImg);
        clear data_sub
        corrImage_tun = getPixelCorrelationImage(data_tun_reg);
        clear data_tun_reg
    else
        corrImage_tun = getPixelCorrelationImage(data_tun_reg);
        clear data_tun_reg
    end
    
    figure;
    suptitle([mouse '-' expDate])
    colormap gray
    subplot 121
    imagesc(corrImage_bx)
    subplot 122
    imagesc(corrImage_tun)
    
    save(fullfile(fnout,'regOuts&Img.mat'),'out_bx','out_tun', 'regImg')
    save(fullfile(fnout,'corrImages.mat'),'corrImage_bx','corrImage_tun')

end

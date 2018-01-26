clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
slct_expt = 4:6 % done: 1,7:8
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
fnout = fullfile(fn,'data processing');
nFrPerRun = [];
input_bx = [];
nFrPerRun = zeros(1,expt(iexp).nrun);
for irun = 1:expt(iexp).nrun

    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
    input = loadMworksFile(SubNum,expDate,expTime);
    nFrPerRun(irun) = size(data_temp,3);
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
            input_bx = [input_bx input];
        end
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);


load(fullfile(fnout,'regOuts&Img.mat'))

out_bx = double(out_bx);
[~,data_bx_reg] = stackRegister_MA(data_bx,[],[],out_bx);
clear data_bx
fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
fName = [dirFolder '_000_000'];

input_tun = loadMworksFile(SubNum,expDate,dirTime);
data_tun = loadsbx_choosepmt(1,mouse,expDate,dirFolder,fName);

% down-sample
data_tun_down = stackGroupProject(data_tun,down);
clear data_tun
nfr_tun = size(data_tun_down,3);

% remove negative data by subtraction
data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
data_tun = data_sub;
clear data_sub

[out_tun,data_tun_reg] = stackRegister(data_tun,regImg);
%%
getTaskMaxDFF
clear data_bx_reg
dirTuningMaxDFF
clear data_tun data_tun_reg
end
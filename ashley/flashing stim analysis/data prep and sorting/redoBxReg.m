clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1

for iexp = [12,13,19,20]

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;

%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun

    if irun == 3 & strcmp(expDate,'150508')
        continue
    end
    
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    [input, data_temp, t] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName);  

    disp(t)
    
    if irun == 1
        data_bx = data_temp;
        input_bx = input;
    else
        data_bx = cat(3,data_bx,data_temp);
        input_bx = [input_bx input];
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx,3);

% remove negative data by subtraction
data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
data_bx = data_sub;
clear data_sub

xpix = size(data_bx,2);
ypix = size(data_bx,1);

fnin = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

load(fullfile(fnin,'regOuts&Img.mat'));

[out_bx data_bx_reg] = stackRegister_MA(data_bx,[],[],out_bx);
clear data_bx


%% ********* behavior max dF/F ************

fnbx = fullfile(rc.ashleyAnalysis,'FSAV Summaries','awFSAVdatasets_V1');
% get behavior max dF/F
taskMaxDFF
end
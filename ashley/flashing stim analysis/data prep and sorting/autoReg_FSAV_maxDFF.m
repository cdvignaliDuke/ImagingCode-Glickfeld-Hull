clear all
close all
%%
ds = '_V1gad';
doDendrites = 0;
rc = behavConstsAV;
eval(['awFSAVdatasets' ds])
slct_expt = 5:size(expt,2);

%%
for iexp = slct_expt
tic
if expt(iexp).greenredsimultaneous == 0
    t(iexp) = toc;
    disp(t)    
    continue
end
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);

%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun
    
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if irun == 3 & strcmp(expDate,'150508')
        [input, data_temp] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName,nframes);  
    else
        [input, data_temp] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName);  
    end
    
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
                    input_bx.(char(genvarname(inpPlus(i)))) = cell(1,80);
                end
            end
            input_temp = [input_bx input];
        end
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);

% remove negative data by subtraction
data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
data_bx = data_sub;
clear data_sub

%% load outs and re-regester data
load(fullfile(fn,'reg2RedOuts&Img.mat'));
[~, data_bx_reg] = stackRegister_MA(data_bx,[],[],double(out_bx));
clear data_bx

%% load tuning data and register
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
down = 10;

fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate, dirFolder);
fName = [dirFolder '_000_000'];  
[input_tun, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);

% down-sample
data_down = stackGroupProject(data,down);
clear data

% remove negative data by subtraction
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% load outs and re-regester data
[out data_tun_reg] = stackRegister_MA(data_sub,[],[],double(out_tun));

%% ********* behavior max dF/F ************

fnbx = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
% get behavior max dF/F
taskMaxDFF

%% ********* tuning max dF/F ************
%%
fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
fnout = fn;
dirTuningMaxDFF



%% save data
t(iexp) = toc;
end
disp(t)
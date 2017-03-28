clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1
for iexp = [7,10,12,13,18,19,20];
tic
SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);

%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun

    if irun == 3 & strcmp(expDate,'150508')
        preRegFrames = size(data_bx,3);
    end
    
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
                    input.(char(genvarname(inpPlus(i)))) = cell(1,80);
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
load(fullfile(fn,'regOuts&Img.mat'));
if strcmp(expDate,'150508')
    [out data_reg1] = stackRegister_MA(data_bx(:,:,1:preRegFrames),[],[],out_bx);
    data_bx = data_bx(:,:,preRegFrames+1:end);
    [out data_reg2] = stackRegister(data_bx,data_corr_img);
    data_reg = cat(3,data_reg1,data_reg2);
    clear data_reg1 data_reg2 data_bx
else
    [out data_reg] = stackRegister_MA(data_bx,[],[],out_bx);
    clear data_bx
end
%% load mask
load(fullfile(fn,'final_mask.mat'))
%% get timecourses
data_tc = stackGetTimeCourses(data_reg,mask_cell);
%% subtract neuropil

% get neuropil timecourses
buf = 4;
np = 6;

data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);
clear data_reg
%% save data
save(fullfile(fn,'timecourses.mat'),'data_tc_subnp')
save(fullfile(fn,'raw_tc.mat'),'out_bx','data_tc','buf','np')
t(iexp) = toc;
disp(t)
end
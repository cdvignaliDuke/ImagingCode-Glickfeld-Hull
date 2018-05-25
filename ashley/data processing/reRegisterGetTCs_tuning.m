clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
doRedChannel = 1;
slct_expt = 2:4;
%%
dataImg = cell(1,length(slct_expt));
dataImgInfo = cell(1,length(slct_expt));
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
fnout = fullfile(fn,dirFolder);
%% load and register all behavior data
% % concatenate data files
% for irun = 1:expt(iexp).nrun
%     runFolder = expt(iexp).dirtuning;
%     expTime = expt(iexp).dirtuning_time;
    fName = [dirFolder '_000_000'];

    data_tun = loadsbx_choosepmt(1,mouse,expDate,dirFolder,fName);
%     if irun == 1
%         data_bx = data_temp;
%     else
%         data_bx = cat(3,data_bx,data_temp);
%     end
%     clear data_temp 
% end
data_tun = stackGroupProject(data_tun,10);
load(fullfile(fnout,'regOuts&Img.mat'))
out_tun = double(out_tun);
[~,data_tun_reg] = stackRegister_MA(data_tun,[],[],out_tun);
dataImg{iexp} = mean(data_tun_reg,3);
dataImgInfo{iexp} = [SubNum '-' expDate];
clear data_tun

%% 
if doRedChannel
    load(fullfile(fnout,'twoColorMask&CellLabels'))
    dataTC = stackGetTimeCourses(data_tun_reg,twoColorMask);
    buf = 4;
    np = 6;
    dataTC_npSub = getWeightedNeuropilTimeCourse(data_tun_reg,dataTC,twoColorMask,buf,np);

    save(fullfile(fnout,'timecourses_tun_cells.mat'),'dataTC_npSub','isLabeledCell')
else
    load(fullfile(fnout,'final_mask.mat'))
    dataTC = stackGetTimeCourses(data_tun_reg,mask_cell);
    buf = 4;
    np = 6;
    dataTC_npSub = getWeightedNeuropilTimeCourse(data_tun_reg,dataTC,twoColorMask,buf,np);
    save(fullfile(fnout,'timecourses_tun_cells.mat'),'dataTC_npSub','dataTC')
end
end
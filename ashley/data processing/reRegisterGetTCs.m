clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
slct_expt = 8%1:size(expt,2);
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
fnout = fullfile(fn,'data processing');
%% load and register all behavior data
% concatenate data files
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

load(fullfile(fnout,'regOuts&Img.mat'))
out_bx = double(out_bx);
[~,data_bx_reg] = stackRegister_MA(data_bx,[],[],out_bx);
dataImg{iexp} = mean(data_bx_reg,3);
dataImgInfo{iexp} = [SubNum '-' expDate];
clear data_bx

%% 
load(fullfile(fnout,'twoColorMask&CellLabels'))
dataTC = stackGetTimeCourses(data_bx_reg,twoColorMask);
buf = 4;
np = 6;
dataTC_npSub = getWeightedNeuropilTimeCourse(data_bx_reg,dataTC,twoColorMask,buf,np);

save(fullfile(fnout,'timecourses_bx_cells.mat'),'dataTC_npSub','isLabeledCell','isFromActiveMask')
end

[nrows,ncols] = optimizeSubplotDim(length(slct_expt));
figure; colormap gray
for iexp = 1:length(slct_expt)
    subplot(nrows,ncols,iexp)
    imagesc(dataImg{iexp})
    title(dataImgInfo{iexp})
end
print(fullfile(rc.caOutputDir, ds, 'avgRegImgEaExpt'),'-dpdf','-fillpage')

figure;colormap gray
for iexp = 1:length(slct_expt)
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
    load(fullfile(fnout,'regOuts&Img.mat'))
    subplot(nrows,ncols,iexp)
    imagesc(regImg)
    title(dataImgInfo{iexp})
end
    
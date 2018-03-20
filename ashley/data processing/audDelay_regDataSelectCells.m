clear all
close all
ds = 'audDelay_V1_EMX';
rc = behavConstsAV;
eval(ds)
slct_expt = 1;
%%
iexp = slct_expt

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;

fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fullfile(fn,'data processing'),'dir')
    mkdir(fn,'data processing')
end
fnout = fullfile(fn,'data processing');
%% load and register all behavior and tuning data
% concatenate data files
nFrPerRun = zeros(1,expt(iexp).nrun);
for irun = 1:expt(iexp).nrun

    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if strcmp(ds, 'audDelay_V1_EMX') & irun == 3 & strcmp(expDate,'180302')
        input = loadMworksFile(SubNum,expDate,expTime,rc.behavData);
        data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName,31174,rc);      
    else
        input = loadMworksFile(SubNum,expDate,expTime,rc.behavData);
        data_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName,[],rc);  
    end
    
    disp(t)
    nFrPerRun(irun) = size(data_temp,3);
    if irun == 1
        data = data_temp;
        mworks = input;
    else
        data = cat(3,data,data_temp);
        try
            mworks = [mworks input];
        catch
            inpNames1 = fieldnames(mworks);
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
                    mworks.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames1));
                end
            end
            mworks = [mworks input];
        end
    end
    clear data_temp input
end
mworks = concatenateDataBlocks(mworks);
nfr_bx = size(data,3);

% remove negative data by subtraction
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

% register 

regImg = mean(data(:,:,expt(iexp).regImgStartFrame),3);
figure;imagesc(regImg);colormap gray
[out_img,data_reg] = stackRegister(data,regImg);
clear data
save(fullfile(fnout,'regOuts&Img.mat'),'out_img', 'regImg')


% save some F images as tifs
nImg = 200;
frInterval = 1:floor((nfr_bx-100)/nImg):nfr_bx-100;
F_bx = nan(size(data_reg,1),size(data_reg,2),nImg);
for iimg = 1:nImg
    F_bx(:,:,iimg) = mean(data_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
end
writetiff(F_bx,fullfile(fnout,'FImages_bx'))
%% get task max dF/F for 1st resp, mid resp, and target resp

maxDFF_ori = getAudDelayMaxDFF(params,mworks,data_reg,nFrPerRun);
%% crop cells image

%crop it
selection_img = max(maxDFF_ori,[],3);
figure;colormap gray; imagesc(selection_img)

[ypix,xpix,~] = size(selection_img);

%**enter vals here***
xcrop = [1:10 790:xpix];
ycrop = [1:10 260:ypix];

% check if behavior image still needs cropping
figure;colormap gray; imagesc(selection_img)
selection_img_crop = selection_img;
selection_img_crop(:,xcrop) = 0;
selection_img_crop(ycrop,:) = 0;

imagesc(selection_img_crop)

selection_img_crop = maxDFF_ori;
selection_img_crop(:,xcrop,:) = 0;
selection_img_crop(ycrop,:,:) = 0;

save(fullfile(fnout,'max_images_crop.mat'),'selection_img_crop','xcrop','ycrop');
close all

%% get cells mask
mask_cell = maskFromMultiMaxDFFStack(selection_img_crop);

figure; setFigParams4Print('portrait')
imagesc(mask_cell);
title({[num2str(length(unique(mask_cell(:)))-1) ' cells'];[mouse '-' expDate]})
print(fullfile(fnout,'final_mask_cells'),'-dpdf')

save(fullfile(fnout,'final_mask_cells.mat'),'mask_cell');
close all
%% get time-courses
buf = 4;
np = 6;

data_tc = stackGetTimeCourses(data_reg,mask_cell);
data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);

save(fullfile(fnout,'timecourses_cells.mat'),'data_tc_subnp','data_tc','buf','np')

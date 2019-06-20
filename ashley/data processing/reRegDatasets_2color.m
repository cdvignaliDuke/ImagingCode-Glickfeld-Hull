clear all
close all
ds = 'FSAV_V1_SOM';
rc = behavConstsAV;
eval(ds)
slct_expt = 13;
doPreviousReg = false;
%%
iexp = slct_expt

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
if strcmp(rc.name, 'ashle')
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end
elseif strcmp(rc.name, 'carolyn')
    fn = fullfile(rc.carolynAnalysis,mouse,'Two-Photon Imaging',expDate);
    if ~exist(fullfile(fn,'data processing'),'dir')
        mkdir(fn,'data processing')
    end  
else
    error('you do not belong here')
end
fnout = fullfile(fn,'data processing');
%% load and register all behavior and tuning data
% concatenate data files
nFrPerRun = zeros(1,expt(iexp).nrun);
for irun = 1:expt(iexp).nrun

    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if ~isempty(expt(iexp).nframesPerRun)
        nframes = expt(iexp).nframesPerRun(irun);
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
        data_temp_r = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName);
%         [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
%             SubNum,expDate,expTime,mouse,runFolder,fName,nframes);
    else
        data_temp_g = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName,expt(iexp).nframesPerRun);
        data_temp_r = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName,expt(iexp).nframesPerRun);
%         [input, data_temp, t] = Load_SBXdataPlusMWorksData(...
%             SubNum,expDate,expTime,mouse,runFolder,fName);
    end
    input = loadMworksFile(SubNum,expDate,expTime);
    
    nFrPerRun(irun) = size(data_temp_r,3);
    if irun == 1
        data_bx_g = data_temp_g;
        data_bx_r = data_temp_r;
        input_bx = input;
    else
        data_bx_g = cat(3,data_bx_r,data_temp_g);
        data_bx_r = cat(3,data_bx_r,data_temp_r);
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
    clear data_temp_g data_temp_r input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx_r,3);

% remove negative data by subtraction
data_sub = data_bx_g-min(min(min(data_bx_r,[],1),[],2),[],3);
data_bx_g = data_sub;
data_sub = data_bx_r-min(min(min(data_bx_r,[],1),[],2),[],3);
data_bx_r = data_sub;
clear data_sub

% register passive expt

% register 

    regImg = mean(data_bx_r(:,:,expt(iexp).regImgStartFrame:(expt(iexp).regImgStartFrame+99)),3);
    figure;imagesc(regImg);colormap gray
    [out_bx,data_bx_r_reg] = stackRegister(data_bx_r,regImg);
    [~,data_bx_g_reg] = stackRegister_MA(data_bx_g,[],[],out_bx);
%     [out_tun,data_tun_reg] = stackRegister(data_tun,regImg);
    clear data_bx_g data_bx_r

% save some F images as tifs
nImg = 200;
frInterval = 1:floor((nfr_bx-100)/nImg):nfr_bx-100;
F_bx = nan(size(data_bx_g_reg,1),size(data_bx_g_reg,2),nImg);
for iimg = 1:nImg
    F_bx(:,:,iimg) = mean(data_bx_g_reg(:,:,frInterval(iimg):frInterval(iimg)+99),3);
end
writetiff(F_bx,fullfile(fnout,'FImages_bx'))
% F_tun = mean(data_tun_reg,3);
%% get task max dF/F for 1st resp, mid resp, and target resp

getTaskMaxDFF

close all
%% crop cells image

%crop it
% tun_img = max(dFF_dirmax,[],3);
bx_img = max(dFF_bxMax,[],3);
% tuning image
% figure;colormap gray; imagesc(tun_img)

figure;colormap gray; imagesc(bx_img)
%**enter vals here***
xcrop = [1:20 780:xpix];
ycrop = [1:10 260:ypix];

% tun_crop = tun_img;
% tun_crop(:,xcrop) = 0;
% tun_crop(ycrop,:) = 0;
% 
% imagesc(tun_crop)

% check if behavior image still needs cropping
bx_crop = bx_img;
bx_crop(:,xcrop) = 0;
bx_crop(ycrop,:) = 0;

imagesc(bx_crop)

% crop orig max images
% dir_crop = dFF_dirmax;
% dir_crop(:,xcrop,:) = 0;
% dir_crop(ycrop,:,:) = 0;

bx_crop = cat(3, start_max, long_max, tar_max);
bx_crop(:,xcrop,:) = 0;
bx_crop(ycrop,:,:) = 0;

save(fullfile(fnout,'max_images_crop.mat'),'bx_crop','xcrop','ycrop');
close all

%% load red image
redImage_bx = mean(data_bx_r_reg,3);
figure; imagesc(redImage_bx);
writetiff(redImage_bx,fullfile(fnout,'redImageBx'))

redFolder = expt(iexp).redChannelRun;
data_red = loadsbx_choosepmt(2,mouse,expDate,redFolder,[redFolder '_000_000']);
[~,data_red_reg] = stackRegister(data_red,regImg);
redImage = mean(data_red_reg,3);
figure; imagesc(redImage)
writetiff(redImage,fullfile(fnout,'redChannelRun'))
writetiff(stackGroupProject(data_red_reg,10),fullfile(fnout,'redChannelStack'))

%% get cells mask

dFF_stack = cat(3,redImage,bx_crop);
[red_mask_cell, green_mask_cell] = maskFromMultiMaxDFFStack_2color(dFF_stack);

figure; setFigParams4Print('portrait')
subplot 211
imagesc(green_mask_cell);
title({sprintf('%s %s- cells with behavior',...
    num2str(length(unique(green_mask_cell(:)))-1), ...
    expt(iexp).redChannelLabel);[mouse '-' expDate]})
subplot 212
imagesc(red_mask_cell);
title({sprintf('%s %s- cells with behavior',...
    num2str(length(unique(red_mask_cell(:)))-1), ...
    expt(iexp).redChannelLabel);[mouse '-' expDate]})
print(fullfile(fnout,'final_mask_cells'),'-dpdf')

save(fullfile(fnout,'final_mask_cells.mat'),'green_mask_cell','red_mask_cell');
close all

%% get time-courses
buf = 4;
np = 6;

data_bx_r_tc = stackGetTimeCourses(data_bx_g_reg,red_mask_cell);
data_bx_r_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_g_reg,data_bx_r_tc,red_mask_cell,buf,np);

data_bx_g_tc = stackGetTimeCourses(data_bx_g_reg,green_mask_cell);
data_bx_g_tc_subnp = getWeightedNeuropilTimeCourse(data_bx_g_reg,data_bx_g_tc,green_mask_cell,buf,np);

if ~exist(fullfile(fn,dirFolder),'dir')
    mkdir(fntun,dirFolder)
end
save(fullfile(fnout,'timecourses_bx_cells.mat'),...
    'data_bx_g_tc_subnp','data_bx_g_tc',...
    'data_bx_r_tc_subnp','data_bx_r_tc','buf','np')

clear data_bx_g_reg data_bx_r_reg
%% register passive expt, get TCs
fName = [expt(iexp).passExpt '_000_000'];
data_pass_g = loadsbx_choosepmt(1,mouse,expDate,expt(iexp).passExpt,fName);
data_pass_r = loadsbx_choosepmt(2,mouse,expDate,expt(iexp).passExpt,fName);
nfr_pass = size(data_pass_g,3);

[out_pass,data_pass_r_reg] = stackRegister(data_pass_r,regImg);
[~,data_pass_g_reg] = stackRegister_MA(data_pass_g,[],[],out_pass);
save(fullfile(fnout,'regOuts&Img.mat'),'out_bx','out_pass','regImg')
clear data_pass_g data_pass_r

frInterval_pass = 1:floor((nfr_pass-100)/nImg):nfr_pass-100;
F_pass = nan(size(data_pass_g_reg,1),size(data_pass_g_reg,2),nImg);
for iimg = 1:nImg
    F_pass(:,:,iimg) = mean(data_pass_g_reg(:,:,frInterval_pass(iimg):frInterval_pass(iimg)+99),3);
end
writetiff(F_pass,fullfile(fnout,'FImages_pass'))


data_pass_r_tc = stackGetTimeCourses(data_pass_g_reg,red_mask_cell);
data_pass_r_tc_subnp = getWeightedNeuropilTimeCourse(data_pass_g_reg,data_pass_r_tc,red_mask_cell,buf,np);

data_pass_g_tc = stackGetTimeCourses(data_pass_g_reg,green_mask_cell);
data_pass_g_tc_subnp = getWeightedNeuropilTimeCourse(data_pass_g_reg,data_pass_g_tc,green_mask_cell,buf,np);

save(fullfile(fnout,'timecourses_pass_cells.mat'),...
    'data_pass_g_tc_subnp','data_pass_g_tc',...
    'data_pass_r_tc_subnp','data_pass_r_tc','buf','np')

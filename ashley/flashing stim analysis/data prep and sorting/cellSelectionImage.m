clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1

%% direction tuning max dF/F needed?
doDirMax = 0;

%%
for iexp = 3:size(expt,2);

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
down = 10;

%%
if doDirMax
    
%% load direction tuning data
expTime = expt(iexp).dirtuning_time;
fName = [dirFolder '_000_000'];

[input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,dirFolder,fName);    

fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,dirFolder);

if ~exist(fnout,'dir')
    mkdir(fnout)
end

%% down sample and register data according to direction tuning data
data_down = stackGroupProject(data,down);
clear data

data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
data = data_sub;
clear data_sub

%*******CHOOSE REG IMG***************
avg_ind_seed = 900;
avg_ind_start = (avg_ind_seed:100:4*avg_ind_seed)';

avg_ind = [avg_ind_start avg_ind_start+10];
figure
colormap gray
for i = 1:4
    subplot(2,2,i)
    data_avg = mean(data(:,:,avg_ind(i,:)),3);
    imagesq(data_avg);
    title_str = num2str(avg_ind(i,:)');
    title([title_str(1,:) '-' title_str(2,:)])
end

%%  select avg frames and register data, or load previously registered data;

%register data
avg_frames = 800:810;
data_reg = regAndSaveData(data,avg_frames,fnout);



%% image specs
xpix = size(data_reg,2);
ypix = size(data_reg,1);
nfr = size(data_reg,3);

%% get max dF/F projection
dirTuningMaxDFF

%% ******choose crop parameters*******
%**enter vals here***
xcrop = [1:2 790:xpix];
ycrop = [1:2 262:ypix];

maxDFF_crop = maxDFF;
maxDFF_crop(:,xcrop) = 0;
maxDFF_crop(ycrop,:) = 0;

imagesc(maxDFF_crop)

%% get max dF/F for each direction
dirTuningMaxDFF_stimSpecific

setFigParams4Print('landscape')
figure;
colormap gray
for istim = 1:nstim
   subplot(2,nstim/2,istim)
   imagesc(dFF_dirmax(:,:,istim))
   title([num2str(dir(istim)) ' deg'])
end

print(fullfile(fnout,'dir_max_images_fig'),'-dpdf','-fillpage')
save(fullfile(fnout,'dir_max_images.mat'),'dFF_dirmax')

end
%% task driven activity images

taskMaxDFF

end
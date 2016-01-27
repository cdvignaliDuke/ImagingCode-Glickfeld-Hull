clear all
awFSAVdatasets_PM
% for iexp = 1:size(expt,2)
    iexp = 2;

% for iexp = 1:size(expt,2);
SubNum = expt(iexp).SubNum;
date = expt(iexp).date;
time = expt(iexp).rettuning{2,:};
ImgFolder = expt(iexp).rettuning{1,:};
mouse = expt(iexp).mouse;
fName = [ImgFolder '_000_000'];

Load_SBXdataPlusMWorksData    
 
down = 10;
nON = input.nScansOn/down;
nOFF = input.nScansOff/down;
nStim = input.gratingElevationStepN.*input.gratingAzimuthStepN;


% data_reg = readtiff('Retinotopy_V1.tif');
data_down = stackGroupProject(data,down);
clear data

data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
data = data_sub;
clear data_sub

% data_avg = mean(data(:,:,2000:2100),3);
% figure; imagesq(data_avg); colormap(gray)

%get direction tuning registration image and get cells
dirFolder = expt(iexp).dirtuning;
fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileDirMasks);
load('regImg.mat');
load('mask&TCDir.mat');
load('neuropil.mat');
clear npTC npSubTC data_TC

[out data_reg] = stackRegister(data, data_avg);
clear data

try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
if (mod(nRep,1)) >0
    nframes = floor(nRep)*((nON+nOFF)*nStim)
    data_reg = data_reg(:,:,1:nframes);
    nRep = size(data_reg,3)./((nON+nOFF)*nStim);
    nTrials = (nStim.*nRep);
end



%%
writetiff(data_reg, 'Retinotopy_reg2dirtuning');

%write tifs for sorted frames
VSR = 1;
run('sortTrialsAvg_writetiffs.m')

%% get timecourses and subtract neuropil
% get timecourses
try
dataTC = stackGetTimeCourses(data_reg, mask_cell);
catch
    dataTC = stackGetTimeCourses(data_reg, mask_boutons);
end

% get neuropil timecourses
buf = 4;
np = 6;
nCells = size(dataTC,2);


npTC = zeros(size(dataTC));
for i = 1:nCells
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

% npTC = stackGetTimeCourses(data_reg,neuropil);


dataTC_mavg = tsmovavg(dataTC,'s',10,1);
npTC_mavg = tsmovavg(npTC,'s',10,1);

% down = 10;
% dataTC_down = reshape(mean(reshape(dataTC',size(dataTC',1),down,size(dataTC',2)/down),2),size(dataTC',1),size(dataTC',2)/down)';
% npTC_down = reshape(mean(reshape(npTC',size(npTC',1),down,size(npTC',2)/down),2),size(npTC',1),size(npTC',2)/down)';

ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(dataTC_mavg-tcRemoveDC(npTC_mavg*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
dataTCsub = dataTC-bsxfun(@times,tcRemoveDC(npTC),np_w);
data_TC = dataTCsub;
clear data_reg

%% vis stim parameters
VSsize = input.gratingDiameterDeg;
VSdirectionDeg = input.gratingDirectionDeg;
tAz = double(cell2mat(input.tGratingAzimuthDeg));
Az = unique(tAz);
[nAz azPos] = histc(tAz,Az);
tEl = double(cell2mat(input.tGratingElevationDeg));
El = unique(tEl);
[nEl elPos] = histc(tEl,El);
if any(tEl == 0) | any(tAz == 0)
    if any(tEl == 0)
        tEl2 = tEl+1;
    else
        tEl2 = tEl;
    end
    if any(tAz == 0)
        tAz2 = tAz+1;
    else
        tAz2 = tAz;
    end
    pos = cart2pol(tAz2,tEl2);
else
    pos = cart2pol(tAz,tEl);
end
[n posN] = histc(pos,unique(pos));
Rets = unique(pos);

AzPos = NaN(1,nStim);
ElPos = NaN(1,nStim);
for i = 1:nStim
    AzPos(i) = unique(tAz(pos == Rets(i)));
    ElPos(i) = unique(tEl(pos == Rets(i)));
end



%% dF/F by trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);
stimON_ind = nOFF+1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end


%% sort data by trial type
dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanRetResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
errbarRets = zeros(size(data_TC,2),nStim);
for i = 1:nStim
    trials = find(pos(:,1:nTrials) == Rets(i));
    dFoverF_meanRetResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
    errbarRets(:,i) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end

figure;
for i = 1:nStim
    plot(dFoverF_meanRetResp(:,10,i));
    hold on
end

%% plot tuning curves
dFoverF_meanOFFRetResp = (squeeze(mean(dFoverF_meanRetResp(1:10,:,:),1)));

RetRespPerCell = (squeeze(mean(dFoverF_meanRetResp(11:end,:,:),1)));

% az x el matrix per cell
resp = squeeze(mean(dFoverFCellsTrials(11:end,:,:),1));
resp_AzElMat = zeros(size(resp,1),nRep,length(El),length(Az));

for iaz = 1:length(Az)
    for iel = 1:length(El)
        resp_AzElMat(:,:,iel,iaz) = resp(:,find(tEl == El(iel) & tAz == Az(iaz)));
    end
end

meanResp_AzElMat = squeeze(mean(resp_AzElMat,2));
errResp_AzElMat = squeeze(std(resp_AzElMat,[],2)/sqrt(double(nRep)));

%%
figure;
imagesc(squeeze(mean(meanResp_AzElMat,1)));
set(gca,'xTick', [1:length(Az)])
set(gca,'xTickLabel', cellfun(@num2str,num2cell(Az),'UniformOutput',false))
set(gca,'yTick', [1:length(El)])
set(gca,'yTickLabel', cellfun(@num2str,num2cell(El),'UniformOutput',false))
xlabel('Az')
ylabel('El')
colorbar

runstr = expt(iexp).runs(1,:);
if expt(iexp).nrun>1
    for irun = 2:expt(iexp).nrun
        runstr = [runstr '-' expt(iexp).runs(irun,:)];
    end
end
var = whos('-file',fullfile('Z:\analysis\',mouse,'two-photon imaging', date, runstr,[mouse '-' date '-' runstr '-comboInputDataTCplusVar.mat']));
load(fullfile('Z:\analysis\',mouse,'two-photon imaging', date, runstr,[mouse '-' date '-' runstr '-comboInputDataTCplusVar.mat']), var(structfind(var,'name','input')).name)
retUsedAz = input.gratingAzimuthDeg;
retUsedEl = input.gratingElevationDeg;

title({'Retinotopy';[SubNum '-' date]; ['(Az,El) = (' num2str(retUsedAz) ',' num2str(retUsedEl) ')'] })

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

% get position used for behavior


print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder,'RetPreferences'), '-dpdf')
%% save tuning info
save(fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder,'Timecourses.mat') ,'data_TC','npTC')
save(fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder,'RetPreferences.mat'),'dFoverFCellsTrials','dFoverF_meanRetResp','dFoverF_meanRetResp','resp_AzElMat','meanResp_AzElMat','errResp_AzElMat')

%%
% end
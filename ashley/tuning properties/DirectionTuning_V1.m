% %load data
% 
% SubNum = '613';
% date = '150511';
% time = '1513';
% ImgFolder = '005';
% mouse = 'AW13';
% % fName = '003_000_000';
% % 
% % % load MWorks file
% % load MWorks file
% CD = ['Z:\data\' mouse '\mworks\' date];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);


% analysis folder
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

% data_reg = readtiff('DirectionTuning_V1.tif');
%used load_SBXdataset_fast.m to load data
%% Parameters
down = 10;
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nStim = double(input.gratingDirectionStepN);

%% downsample and register data

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,100:110),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

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
%save data_reg
writetiff(data_reg, 'DirectionTuning_V1');


%%
%write tifs for sorted frames
VSR = 2;
run('sortTrialsAvg_writetiffs.m')
%% create dF/F stack

nOFF_ind = zeros(1,(nOFF*nStim*nRep));
nOFF_1 = zeros(1,nRep*nStim);
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    nOFF_1(1,iStim) = 1+((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

%% movie with trial start marked
data_offmark = data_reg;
data_offmark(1:55,721:end,nOFF_1) = max(max(mean(data_reg,3)));
writetiff(data_offmark,'rawFtrialstartmark')

%%

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

% dF average F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF = bsxfun(@rdivide,dF_data,nOFF_avg);

%motion/global F change index
motionThreshold = 0.05;
ind_motion = find(abs(diff(squeeze(mean(mean(dFoverF,1),2)),1))>motionThreshold);

% writetiff(dFoverF(:,:,ind_motion+1),'motionFrames');

m = setdiff(1:size(data_reg,3),ind_motion+1);

% maxF = max(data_reg(:,:,nON_ind),[],3);
% meanF_on = mean(data_reg(:,:,nON_ind),3);
% meanF_off = mean(data_reg(:,:,nOFF_ind),3);
% figure;imagesq(maxF); colormap(gray)
% figure;imagesq(meanF_on); colormap(gray);title('F on')
% figure;imagesq(meanF_off);colormap(gray);title('F off')
% sub = bsxfun(@minus,meanF_on,meanF_off);
% sub(sub < 0) = 0;
% figure;imagesq(sub);colormap(gray);title('on-off')
max_dF = max(dF_data,[],3);
maxDFoverF = max(dFoverF(:,:,m),[],3);
meanDFoverF = mean(dFoverF,3);
maxDFoverFcrop = maxDFoverF;
maxDFoverFcrop(:,[1:32 758:end]) =0;
figure; imagesq(maxDFoverF); colormap(gray)
% figure; imagesq(meanDFoverF); colormap(gray)

%get rid of bright patch to make gui easier to use
% maxDFoverF(:,750:end,:) = 0;

%save max DF/F
writetiff(maxDFoverF, 'maxDFoverF');

bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)

% xcalib = size(data_reg,2)/500;
% ycalib = size(data_reg,1)/250;
% scalebar50umx = 50*xcalib;
% 
% scalebarimg = zeros(size(max_dF));
% scalebarimg(size(data_reg,1)-60-5:size(data_reg,1)-60,size(data_reg,2)-200-scalebar50umx:size(data_reg,2)-200) = 0.5;
% figure;imagesq(scalebarimg);colormap(gray)
% writetiff(scalebarimg,'scalebar')
% 
% max_dF_wSB = cat(3,max_dF,scalebarimg);
% writetiff(max_dF,'maxDF')

%% specific cell mask
% mask_cell33 = double(mask_cell == 33);
% mask_cell121 = double(mask_cell == 121);
% figure;imagesq(mask_cell121);colormap(gray)
% 
% writetiff(mask_cell33,'mask_cell33')
% writetiff(mask_cell121,'mask_cell121')

%% use correlation dF/F to find ROIS

% %use max dF if too many cells
% % b = 5;
% % siz = size(data_reg);
% % corr_map = zeros(siz(1),siz(2));
% % for ix = b:siz(2)-b
% %     for iy = b:siz(1)-b
% %         TC = data_reg(iy,ix,:);
% %         surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
% %         R = corrcoef(TC,surround);
% %         corr_map(iy,ix) = R(1,2);
% %     end
% % end
% % figure; imagesq(corr_map); colormap(gray)
% 
% bwout = imCellEditInteractive(corr_map);
% mask_cell = bwlabel(bwout);

% %timecourses
% data_TC = stackGetTimeCourses(data_reg,mask_cell);
% figure; tcOffsetPlot(data_TC)
%%
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('mask&TCDir.mat','mask_cell','data_TC');
save('regImg.mat', 'data_avg');
% load('mask&TCDir.mat')
%% neurpil subtraction
nCells = size(data_TC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);

npTC = zeros(size(data_TC));
for i = 1:size(data_TC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
npSubTC = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('neuropil.mat','neuropil','npTC','npSubTC');

data_TC = npSubTC;
%%
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;
nRep = size(data_TC,1)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);

%% dF/F by trial
% F per trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end
%% dF/F (by cell) for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i));
    dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp(:,1,i));
    hold on
end

%% find magnitude of response to stim
dFoverF_meanOFFDirResp = (squeeze(mean(dFoverF_meanDirResp(1:10,:,:),1)));
DirRespPerCell = (squeeze(mean(dFoverF_meanDirResp(11:end,:,:),1)));

figure;
plot(DirRespPerCell(1,:))

%% find direction preference
DirRespPerCell_sub = DirRespPerCell-min(min(min(DirRespPerCell,[],1),[],2));
[dirPref_val,dirPref_ind] = max(DirRespPerCell_sub,[],2);

pref2orth_dir = [nStim/2+1:nStim 1:nStim/2];
dirOrth_ind = zeros(size(dirPref_ind));
for i = 1:nStim
    cells = find(dirPref_ind == i);
    dirOrth_ind(cells) = pref2orth_dir(i);
end

dirOrth_val = zeros(size(dirPref_ind)); 
for i = 1:size(data_TC,2)
    ind = dirOrth_ind(i);
    cell = DirRespPerCell_sub(i,:);
    dirOrth_val(i) = cell(ind);
end

cellDSI = zeros(size(dirPref_ind));
cellDSI = (dirPref_val - dirOrth_val)./(dirPref_val + dirOrth_val);
dirSlctvCells = find(cellDSI > 0.3);
figure;hist(cellDSI,10);title('DSI')

%% find orientation preference
dFoverF_meanOriResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim/2);
for i = 1:nStim/2
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i) | DirectionDeg(:,1:nTrials) == Dirs(i+nStim/2));
    dFoverF_meanOriResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end
dFoverF_meanOFFOriResp = (squeeze(mean(dFoverF_meanOriResp(1:10,:,:),1)))+1;
dFoverF_meanONOriResp = (squeeze(mean(dFoverF_meanOriResp(11:end,:,:),1)))+1;
OriRespPerCell = dFoverF_meanONOriResp;

figure;
plot(OriRespPerCell(1,:))

OriRespPerCell_sub = OriRespPerCell-min(min(min(OriRespPerCell,[],1),[],2));
[oriPref_val,oriPref_ind] = max(OriRespPerCell_sub,[],2);

pref2orth_ori = [(nStim/2)/2+1:nStim/2 1:(nStim/2)/2];
oriOrth_ind = zeros(size(oriPref_ind));
for i = 1:nStim/2
    cells = find(oriPref_ind == i);
    oriOrth_ind(cells) = pref2orth_ori(i);
end

oriOrth_val = zeros(size(oriPref_ind)); 
for i = 1:size(data_TC,2)
    ind = oriOrth_ind(i);
    cell = OriRespPerCell_sub(i,:);
    oriOrth_val(i) = cell(ind);
end

cellOSI = zeros(size(oriPref_ind));
cellOSI = (oriPref_val - oriOrth_val)./(oriPref_val + oriOrth_val);
oriSlctvCells = find(cellOSI > 0.3);
figure;hist(cellOSI,10);title('OSI')

%% get tuning curves

dFoverFDirResp = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp(i,:) = squeeze(mean(mean(dFoverFCellsTrials(11:16,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end

%% save tuning info
save('TuningPreferences.mat','oriPref_ind','dirPref_ind','dirSlctvCells','oriSlctvCells','dFoverFDirResp','dFoverF_meanDirResp')

%% plot cell tuning
% plot responses to each direction
cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cMap = colormap(jet(nStim));
start = 1;
errbar = zeros(nStim,size(data_TC,2));
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:31
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(8,4,iplot);
    for i = 1:nStim
        plot(dFoverF_meanDirResp(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['Cell ' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp(:,cell,i),[],1);
        if ymax_i > ymax
            ymax = ymax_i;
        end
        axis([0 20 -0.05 ymax]);
        hold on
        legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
        hold on

    end   
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(8,4,iplot),'color',[0.9 0.9 0.9])
    end
    if ismember(cell,cellsPrefNinety) > 0
        set(subplot(8,4,iplot),'color',[0.8 0.8 0.8])
    end
    if iplot == 15 & ifig == 1
        legend(legendInfo,'Location','SouthEast')
    end
    start = start+1;
    hold on
end
end

% plot tuning curves
start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:31
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(8,4,iplot);
    errorbar(dFoverFDirResp(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp(:,cell),[],1);
    if ymax_i > ymax
        ymax = ymax_i;
    end
    title(['Cell ' num2str(cell)]);
    hold on
    axis([0 length(Dirs) -0.05 ymax]);
    hold on
%     if ismember(cell,baselineStimRespIndex_V) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
    hold on
    start = start+1;
end
end




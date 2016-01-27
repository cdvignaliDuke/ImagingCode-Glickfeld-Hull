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
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;

%% downsample and register data

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,300:310),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%save data_reg
writetiff(data_reg, 'DirectionTuning_axons');

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
%%
%write tifs for sorted frames
VSR = 2;
run('sortTrialsAvg_writetiffs.m')

%% Parameters
pre_win = [floor(0.75*nOFF) nOFF];
post_win = [nOFF+1 nOFF+nON];

dirs = unique(cell2mat(input.tGratingDirectionDeg));
nCond = length(dirs);

alpha = .05./(nCond);
b= 5;
pix = 5;
bouton_diam = 5;

%% sort data
nRep_mat = zeros(1,nCond);
siz = size(data_reg);
data_sort = zeros(siz);
start = 1;
for iCond = 1:nCond
    ids = find(cell2mat(input.tGratingDirectionDeg) == dirs(iCond));
     nRep_mat(1,iCond) = length(ids);
    for iRep = 1:length(ids)
        trial = ids(iRep);
        data_sort(:,:,start:start+nON+nOFF-1) = data_reg(:,:,1+((trial-1)*(nON+nOFF)):trial*(nON+nOFF));
        start = start+nON+nOFF;
    end
end

%% create resp structure
resp = struct;
start = 1;
for iCond = 1:nCond; 
    resp(iCond).on = [];
    resp(iCond).off = [];
    roi_stim = zeros(siz(1),siz(2),nOFF+nON,nRep_mat(1,iCond));
    for iRep = 1:nRep_mat(1,iCond);
        roi_stim(:,:,:,iRep) = data_sort(:,:,start:start-1+nON+nOFF);
        start = start+nON+nOFF;
    end
    resp(iCond).off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
    resp(iCond).on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));
end

%% pixel based ttest
Info_ttest_mat = zeros(siz(1),siz(2),nCond);

f1 = fspecial('average');
for iCond = 1:nCond
    resp(iCond).on_long = reshape(resp(iCond).on, [siz(1) siz(2)*nRep_mat(1,iCond)]);
    resp(iCond).off_long = reshape(resp(iCond).off, [siz(1) siz(2)*nRep_mat(1,iCond)]);
    resp(iCond).on_sm = reshape(filter2(f1,resp(iCond).on_long),[siz(1) siz(2) nRep_mat(1,iCond)]);
    resp(iCond).off_sm = reshape(filter2(f1,resp(iCond).off_long),[siz(1) siz(2) nRep_mat(1,iCond)]);
end

for iy = b+1:siz(1)-b
    fprintf([num2str(iy) ' '])
    for ix = b+1:siz(2)-b
        p_ttestB = zeros(1,1,nCond);
        for iCond = 1:nCond;
            [h_ttestB1,p_ttestB1] = ttest(resp(iCond).off_sm(iy,ix,:),resp(iCond).on_sm(iy,ix,:),alpha,'left');
            p_ttestB(1,1,iCond) = p_ttestB1;
        end
    Info_ttest_mat(iy,ix,:) = p_ttestB;
    end
end

Info_ttest_mat_long = interp2(reshape(Info_ttest_mat, [siz(1) siz(2)*nCond]));
Info_ttest_smooth = filter2(f1,Info_ttest_mat_long);
siz_interp = size(Info_ttest_smooth);
Xi = 1:2:siz_interp(1);
Yi = 1:2:siz_interp(2);
ttest_smooth_siz = interp2(Info_ttest_smooth, Yi', Xi);
ttest_smooth = min(reshape(ttest_smooth_siz, [siz(1) siz(2) nCond]),[],3);
ttest_mask = min(ttest_smooth,[],3) < alpha;

ttest_mask(1:b,1:end) = 0;
ttest_mask(1:end, 1:b) = 0;
ttest_mask(1:end, siz(2)-b:end) = 0;
ttest_mask(siz(1)-b:end,1:end) = 0;

figure; imagesq(ttest_mask);
        
%% find dF/F for iCond to find active boutons
resp_dFoverF = zeros(siz(1),siz(2),nCond);
for iCond = 1:nCond
    resp_dFoverF(:,:,iCond) = (mean(resp(iCond).on,3)-mean(resp(iCond).off,3))./(mean(resp(iCond).off,3));
end

max_dF = max(resp_dFoverF,[],3);
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS- local maxima

max_interp = interp2(max_dF);
f1 = fspecial('average');   
max_interp_sm = filter2(f1, max_interp);
siz2 = size(max_interp);
Xi = 1:2:siz2(1);
Yi = 1:2:siz2(2);
stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);

local_max = zeros(siz(1), siz(2));
for iy = b+1:(siz(1)-b);
    for ix = b+1:(siz(2)-b);            
        sub = stack_max_interp_sm(iy-pix:iy+pix,ix-pix:ix+pix);
        sub_long = reshape(sub, [1 (pix*2+1)^2]);
        [sub_long_order ind_order] = sort(sub_long);
        if ind_order(end)==ceil(((pix*2+1)^2)/2)
            local_max(iy,ix) = 1;
        end
    end
end
ind_local_max = find(local_max == 1);
figure; imagesq(local_max);

%% combine ttest and local maxima
ttest_long = reshape(ttest_smooth, [siz(1)*siz(2) 1]);
ind_highP = find(ttest_long(ind_local_max,:)>=(alpha));
local_max_long = reshape(local_max, [siz(1)*siz(2) 1]);
local_max_long(ind_local_max(ind_highP,:),:) = 0;
local_max_sig = reshape(local_max_long, [siz(1) siz(2)]);

n_pix = sum(sum(local_max_sig));
[i, j] = find(local_max_sig ==1);

%expand local maxima
pix_add = floor(bouton_diam/2);
FOV = zeros(size(local_max_sig));
for ipix = 1:n_pix
    FOV(i(ipix)-pix_add:i(ipix)+pix_add,j(ipix)-pix_add:j(ipix)+pix_add) = 1;
end

figure; imagesq(FOV)

%% Get timecourses

mask_boutons = bwlabel(FOV);
figure; imagesq(mask_boutons);
data_TC = stackGetTimeCourses(data_reg,mask_boutons);

figure; tcOffsetPlot(data_TC)

%% save raw timecourse and mask
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('mask&TCDir.mat','mask_boutons','data_TC');
save('regImg.mat', 'data_avg');

%% neurpil subtraction
nBoutons = size(data_TC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_boutons,buf,np);

npTC = zeros(size(data_TC));
for i = 1:size(data_TC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nBoutons);
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

%% find direction and orientation tuning of boutons
% additonal parameters:

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
    plot(dFoverF_meanDirResp(:,4,i));
    hold on
end

%% find magnitude of response to stim

dFoverF_meanOFFDirResp = (squeeze(mean(dFoverF_meanDirResp(1:10,:,:),1)));

dFoverF_meanONDirResp = (squeeze(mean(dFoverF_meanDirResp(11:end,:,:),1)));
DirRespPerCell = dFoverF_meanONDirResp;
% DirRespPerCell = dFoverF_meanONDirResp./dFoverF_meanOFFDirResp;

figure;
plot(DirRespPerCell(4,:))

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
figure;hist(cellDSI,10)

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
plot(OriRespPerCell(4,:))

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
figure;hist(cellOSI,10);

%% get tuning curves

dFoverFDirResp = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp(i,:) = squeeze(mean(mean(dFoverFCellsTrials(11:16,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end

%% save tuning info
save('TuningPreferences.mat','oriPref_ind','dirPref_ind','dirSlctvCells','oriSlctvCells','dFoverFDirResp')

%% plot cell tuning
% % plot responses to each direction
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
% cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
% cMap = colormap(jet(nStim));
% start = 1;
% errbar = zeros(nStim,size(data_TC,2));
% for ifig = 1:ceil(size(data_TC,2)/15)
% figure;
% for iplot = 1:31
% %     cell = baselineStimRespIndex_V(start);
%     cell = start;
%     ymax = .1;
%     subplot(8,4,iplot);
%     for i = 1:nStim
%         plot(dFoverF_meanDirResp(:,cell,i),'color',cMap(i,:));
%         hold on
%         vline(10,'k');
%         hold on
%         title(['Cell ' num2str(cell)]);
%         hold on
%         ymax_i = max(dFoverF_meanDirResp(:,cell,i),[],1);
%         if ymax_i > ymax
%             ymax = ymax_i;
%         end
%         axis([0 20 -0.05 ymax]);
%         hold on
%         legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
%         hold on
% 
%     end   
%     if ismember(cell,cellsPrefZero) > 0
%         set(subplot(8,4,iplot),'color',[0.9 0.9 0.9])
%     end
%     if ismember(cell,cellsPrefNinety) > 0
%         set(subplot(8,4,iplot),'color',[0.8 0.8 0.8])
%     end
%     if iplot == 15 & ifig == 1
%         legend(legendInfo,'Location','SouthEast')
%     end
%     start = start+1;
%     hold on
% end
% end

% % plot tuning curves
% start = 1;
% for ifig = 1:ceil(size(data_TC,2)/15)
% figure;
% for iplot = 1:31
% %     cell = baselineStimRespIndex_V(start);
%     cell = start;
%     ymax = .1;
%     subplot(8,4,iplot);
%     errorbar(dFoverFDirResp(:,cell),errbar(:,cell),'k');
%     hold on
%     ymax_i = max(dFoverFDirResp(:,cell),[],1);
%     if ymax_i > ymax
%         ymax = ymax_i;
%     end
%     title(['Cell ' num2str(cell)]);
%     hold on
%     axis([0 length(Dirs) -0.05 ymax]);
%     hold on
% %     if ismember(cell,baselineStimRespIndex_V) > 0
% %         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
% %     end
%     hold on
%     start = start+1;
% end
% end


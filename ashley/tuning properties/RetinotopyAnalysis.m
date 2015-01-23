SubNum = '607';
date = '150121';
time = '1832';
ImgFolder = '001';
mouse = 'AW07';
fName = '001_000_000';



% data_reg = readtiff('Retinotopy_V1.tif');

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data

% register
data_avg = mean(data_sub(:,:,900:910),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%save data_reg
writetiff(data_reg, 'Retinotopy_V1');

%registered image
data_avg = mean(data_reg(:,:,:),3);
figure; imagesq(data_avg); colormap(gray)

data_reg = readtiff('Retinotopy_V1.tif');
data_reg = data_reg(:,:,1:5400);
% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);


%%

orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = input.nScansOn;
nOFF = input.nScansOff;
nStim = input.gratingElevationStepN.*input.gratingAzimuthStepN;

nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);

%% create dF/F stack
%find off and on frames
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

%dF/F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, nOFF_avg);
% % max_dF = max(dFoverF_data,[],3);
% max_dF = max(dF_data,[],3);
% figure; imagesq(max_dF); colormap(gray)

%find cells with correlation matrix
b = 5;
siz = size(data_reg);
corr_map = zeros(siz(1),siz(2));
for ix = b:siz(2)-b
    for iy = b:siz(1)-b
        TC = data_reg(iy,ix,:);
        surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
        R = corrcoef(TC,surround);
        corr_map(iy,ix) = R(1,2);
    end
end

figure; imagesq(corr_map); colormap(gray)

%% use max dF/F to find ROIS

% bwout = imCellEditInteractive(max_dF);
% mask_cell = bwlabel(bwout);
bwout = imCellEditInteractive(corr_map);
mask_cell = bwlabel(bwout);
    %found with 0.9 and 3pixels

%timecourses
data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

%save
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('mask&TCRet.mat', 'data_TC', 'mask_cell');
%%
% find on indices for the first frame of each stimulus start period and iti (Off) period
for itrial = 1:(nStim*nRep)
    nON_ind_firsts(itrial) = nON_ind(1+(nON*(itrial-1)));
end
for itrial = 1:(nStim*nRep)
    nOFF_ind_firsts(itrial) = nOFF_ind(1+(nOFF*(itrial-1)));
end

%average like trials
nCells = size(data_TC,2);
data_TCtrials = zeros(nON+nOFF,nCells,nTrials);

start = 1;
for itrial = 1:nTrials
    data_TCtrials(1:nON+nOFF,:,itrial) = data_TC(start:start+nOFF+nON-1,:);
    start = start+nOFF+nON;
end

data_TCtrialsSort = zeros(nON+nOFF,nCells,nRep,nStim);
for istim = 1:nStim
    for irep = 1:nRep
        data_TCtrialsSort(:,:,irep,istim) = data_TCtrials(:,:,((irep-1)*nRep)+istim);
        x = ((irep-1)*nRep)+istim
    end
end

data_TCavgTrialType = squeeze(mean(data_TCtrialsSort,3));
for istim = 1:nStim
figure;
for icell = 1:16
    subplot(4,4,icell);
    hold on
        plot(data_TCavgTrialType(nOFF-5:end,icell,istim));
        hold on
    axis([0 (nOFF-(nOFF-5)+nON) -0.1 0.2]);
    hold on
    vline(nOFF-(nOFF-5),'c');
    hold on
end
end

data_TCavgTrialsCells = squeeze(mean(mean(data_TCtrialsSort(:,:,:,:),3),2));
figure;
for istim = 1:nStim
    subplot(2,3,istim)
    hold on
    plot(data_TCavgTrialsCells(:,istim))
    axis([0 (nOFF+nON) -0.05 0.05]);
    hold on
    vline(nOFF,'c');
    hold on
end

%average each trial response per cell
nCells = size(data_TC,2);
dFoverF_meanONbyCell = zeros(nTrials,nCells);
for icell = 1:nCells
    for itrial = 1:nTrials
        tri = data_TC((nON_ind_firsts(itrial):(nON_ind_firsts(itrial)+nON)-1),icell);
        dFoverF_meanONbyCell(itrial,icell) =  mean(tri,1);
    end
end

% resp trace for each trial 
dFoverF_RspbyCellbyTrial = zeros(nON+nOFF,nCells,nTrials);
for icell = 1:nCells
    for itrial = 1:nTrials
        tri = data_TC((nOFF_ind_firsts(itrial):(nON_ind_firsts(itrial)+nON)-1),icell);
        dFoverF_RspbyCellbyTrial(:,icell,itrial) = tri;
    end
end

%remove negative values


%% 

% resp matrix for each position
retResp_mat = zeros(nON+nOFF,nCells,nStim,nRep);

for icell = 1:nCells
    start = 1;
    for istim = 1:nStim
        for irep = 1:nRep
            ind = (nStim*(irep-1)+1)+(istim-1);
            retResp_mat(:,icell,istim,irep) = dFoverF_RspbyCellbyTrial(:,icell,ind);
        end
    end
end

%plot
figure;
for istim = 1:nStim
	subplot(input.gratingAzimuthStepN,input.gratingElevationStepN,istim);
    for icell = 1:nCells
        plot(mean(retResp_mat(:,icell,istim,:),4));
        hold on
    end
    vline(13);
    hold on
end


% figure;
% start = 1;
% for icell = 1:5
%     for istim = 1:nStim
%         pos = istim+(nStim*(start-1));
%         subplot(5,6,pos)
%         for irep = 1:nRep
%             plot(retResp_mat(:,icell,istim,irep));
%             hold on
%         end
%     end
%     start = start+1;
% end
% 
% figure;
% col_mat = strvcat('k', 'b', 'c', 'r', 'm', 'g');
% for icell = 1:nCells
%     subplot(5,6,icell);
%     for istim = 1:nStim
%         plot(retResp_mat(:,icell,istim,1),col_mat(istim,:));
%         hold on
%     end
% end


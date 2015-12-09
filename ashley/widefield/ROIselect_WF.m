%% creates ROIs from retinotopy experiment- only needs to be done once for each mouse
%run iXXX_paths.m (for your mouse) before running code. 
%will also need to manually choose number of ROIs (is hard coded below) and
%assign them to areas 'area_list'

%% Read files and load image stack
cd(fullfile(data_folder, [roi_date '_' mouse], [mouse roi_run]));
imageStack = [mouse roi_run '_MMStack.ome'];
roi_data = double(readtiff([imageStack '.tif']));
roi_data_avg = mean(roi_data,3);
writetiff(roi_data_avg, fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_data_avg.tif']));
cd(behav_folder);
load(['data-i' mouse '-' roi_date '-' roi_time]);

%% Calculate DF/F for each direction in one cycle based on mean of the offs for that direction

nOns = input.nScansOn;
nOffs = input.nScansOff;
nStim = input.gratingDirectionStepN;
nRep = size(input.counterValues,2)./nStim;
nTrials = size(input.counterValues,2);

nTotalFrames = (nOns + nOffs)*nStim*nRep;
siz = size(roi_data);

%make a single stim average movie
trial_data = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_dF = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_dFoverF = zeros(siz(1), siz(2), nOns+nOffs, nTrials);
trial_F = zeros(siz(1), siz(2), nTrials);
for i = 1:nTrials
    trial_data(:,:,:,i) = roi_data(:,:,1+((i-1)*(nOns+nOffs)):i*(nOns+nOffs));
    trial_F(:,:,i) = mean(trial_data(:,:,1+nOffs/2:nOffs,i),3);
    trial_dF(:,:,:,i) = bsxfun(@minus, trial_data(:,:,:,i), trial_F(:,:,i));
    trial_dFoverF(:,:,:,i) = bsxfun(@rdivide, trial_dF(:,:,:,i), trial_F(:,:,i));
end
avg_dFoverF = mean(trial_dFoverF,4);
avg_resp = mean(avg_dFoverF(:,:,nOffs+2:nOffs+10),3);
%% Convert vessel traces on imageJ to vessel map and subtract from ROI mask

vasc_map = readtiff(fullfile(anal_pn,mouse, [mouse '_' roi_date '_vessel_trace.tif']));
figure; imagesc(vasc_map);
vasc_mask = zeros(size(vasc_map)); vasc_mask(find(vasc_map > 60)) = 1;
figure; imagesc(vasc_mask); colormap(gray)

%% Subtract vasculature from DF image, select and store ROI
figure; imagesc(reg_norm); colormap(gray);
clim([0 .9])
%adjust number of ROIs to choose
nROI = 8;
for i = 1:nROI
    roi(i) = impoly;
end

siz = size(avg_resp);
mask = zeros(siz(1),siz(2),nROI);
for i = 1:nROI
    mask(:,:,i) = createMask(roi(i));
end

roi_cluster = sum(mask,3);
mask_cell = bwlabel(roi_cluster);
figure; imagesc(mask_cell);

%adjust list of areas to track
area_list = strvcat('LM','V1','VC','AL','RS','RL','PM','S1');
ind = find(vasc_mask); 
mask_cell_V = mask_cell;
mask_cell_V(ind) = 0;
figure; imagesc(mask_cell_V);

%neuropil mask
np = imCellNeuropil(mask_cell, 2, 5);
np_V = zeros(size(np));
for i = 1:nROI
    np_V_temp = np(:,:,i);
    np_V_temp(ind) = 0;
    np_V(:,:,i) = np_V_temp;
end
figure; imagesc(max(np_V,[],3))
save(fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_masks.mat']), 'roi_cluster', 'mask_cell', 'area_list', 'mask_cell_V', 'np_V');

%% get timecourses from all ROIS with retinotopy
roiTC = stackGetTimeCourses(roi_data, mask_cell_V);
nROI = size(roiTC,2);
roiTC_np = zeros(size(roiTC));
for i = 1:nROI
    roiTC_np(:,i) = stackGetTimeCourses(roi_data, np_V(:,:,i));
end
%roiTC_npsub = roiTC-roiTC_np;


roiTC_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTC_F_stim = zeros(1, nROI, nRep, nStim);
roiTC_DF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTC_DFoverF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_F_stim = zeros(1, nROI, nRep, nStim);
roiTCnp_DF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
roiTCnp_DFoverF_stim = zeros(nOns+nOffs, nROI, nRep, nStim);
for i = 1:nStim
    for ii = 1:nRep
        roiTC_stim(:,:,ii,i) = roiTC(1+((i-1)*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))):(i*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))),:);
        roiTC_F_stim(:,:,ii,i) = mean(roiTC_stim((nOffs/2):nOffs,:,ii,i));
        roiTC_DF_stim(:,:,ii,i) = bsxfun(@minus, roiTC_stim(:,:,ii,i), roiTC_F_stim(:,:,ii,i));
        roiTC_DFoverF_stim(:,:,ii,i) = bsxfun(@rdivide, roiTC_DF_stim(:,:,ii,i), roiTC_F_stim(:,:,ii,i));
        roiTCnp_stim(:,:,ii,i) = roiTC_np(1+((i-1)*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))):(i*(nOns+nOffs))+((ii-1)*(nStim*(nOns+nOffs))),:);
        roiTCnp_F_stim(:,:,ii,i) = mean(roiTCnp_stim((nOffs/2):nOffs,:,ii,i));
        roiTCnp_DF_stim(:,:,ii,i) = bsxfun(@minus, roiTCnp_stim(:,:,ii,i), roiTCnp_F_stim(:,:,ii,i));
        roiTCnp_DFoverF_stim(:,:,ii,i) = bsxfun(@rdivide, roiTCnp_DF_stim(:,:,ii,i), roiTCnp_F_stim(:,:,ii,i));
    end
end
roiTC_DFoverF_avg = squeeze(mean(reshape(roiTC_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTC_DFoverF_sem = squeeze(std(reshape(roiTC_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));

roiTC_DF_avg = squeeze(mean(reshape(roiTC_DF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTC_DF_sem = squeeze(std(reshape(roiTC_DF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));

tt= (1-nOffs:nOns)./(roi_rate./1000);

roiTCnpsub_DFoverF_stim = roiTC_DFoverF_stim-roiTCnp_DFoverF_stim;
roiTCnpsub_DFoverF_avg = squeeze(mean(reshape(roiTCnpsub_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),3));
roiTCnpsub_DFoverF_sem = squeeze(std(reshape(roiTCnpsub_DFoverF_stim, [nOns+nOffs, nROI, nRep*nStim]),[],3)./sqrt(double(nRep*nStim)));

figure;
for i = 1:nROI
    subplot(3,3,i)
    shadedErrorBar(tt, roiTCnpsub_DFoverF_avg(:,i), roiTCnpsub_DFoverF_sem(:,i));
    ylim([-0.01 0.02])
    hold on
    vline(0)
    title(area_list(i,:))
    ylabel('dF/F')
end
subplot(3,3,nROI+1)
plot(tt, roiTCnpsub_DFoverF_avg)
ylim([-.01 0.02])
ylabel('dF/F')
hold on
vline(0)
%legend(area_list, 'Location', 'BestOutside')
print(fullfile(anal_pn, mouse, [mouse '_' roi_date '_allAreaRetResp.pdf']),'-dpdf')

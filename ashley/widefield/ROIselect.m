%% Read files and load image stack
cd('S:\Data\633\Widefield imaging\151103_633\633roi_4directions_33exp_10Hz_1');
imageStack = '633roi_4directions_33exp_10Hz_1_MMStack.ome';
roi_data = double(readtiff([imageStack '.tif']));

mouse = '633';
date = '151103';
time = '1346';
cd('Z:\home\andrew\Behavior\Data');
load(['data-i' mouse '-' date '-' time]);

%% Calculate DF/F for each direction in one cycle based on mean of the offs for that direction

nOns = 50;
nOffs = 50;
nStim = 4;
nRep = 10;
nTotalFrames = (nOns + nOffs)*nStim*nRep;
siz = size(roi_data);
offs_data = cell(1,nStim);
ons_data = cell(1,nStim);
oneD_dataF = zeros(siz(1),siz(2),nOns,nStim);
dataDF = cell(1,nStim);
dataDFoverF = cell(1,nStim);
cycle_avg_DFoverF = zeros(siz(1),siz(2),nStim);
cycle_avg_DF = zeros(siz(1),siz(2),nStim);

k = 1;
for i = 1:nStim
    offs_data{i} = roi_data(:,:,(k*nOffs) - 49:(k*nOffs));
    oneD_dataF(:,:,i) = mean(roi_data(:,:,(k*nOffs) - 49:(k*nOffs)),3);
    k = k + 2;
    ons_data{i} = roi_data(:,:,(i*2*nOns) - 49:(i*2*nOns));    
end

for i = 1:nStim
    dataDF{i} = bsxfun(@minus, ons_data{i}, oneD_dataF(:,:,i));
end

for i = 1:nStim
    dataDFoverF{i} = bsxfun(@rdivide, dataDF{i}, oneD_dataF(:,:,i));
    cycle_avg_DFoverF(:,:,i) = mean(dataDFoverF{i}(:,:,:),3);
    cycle_avg_DF(:,:,i) = mean(dataDF{i}(:,:,:),3);
end

%writetiff(cycle_avg_DFoverF, 'S:\Analysis\633\Widefield imaging\151103\ROI.tif');

%% Convert vessel traces on imageJ to vessel map and subtract from ROI mask

vasc_map = readtiff('C:\Users\shiva\Desktop\vessel_trace.tif');
figure; imagesc(vasc_map);
vasc_mask = zeros(size(vasc_map)); vasc_mask(find(vasc_map > 60)) = 1;
figure; imagesc(vasc_mask); colormap(gray)

%% Subtract vasculature from DF image, select and store ROI

max_DFoverF = max(cycle_avg_DFoverF,[],3);
max_DF = max(cycle_avg_DF,[],3);
%siz_DF = size(max_DF);
%nrow = siz_DF(1);
%ncol = siz_DF(2);
%for i = 1:ncol
%    max_DF_noVasc(:,i) = max_DF(:,i) - vasc_mask(:,i);
%end

figure; imagesc(max_DF); colormap(gray);
for i = 1:6
    roi(i) = impoly;
end

maskV = createMask(roi(1));
maskLM = createMask(roi(2));
maskAL = createMask(roi(3));
maskRL = createMask(roi(4));
maskCV = createMask(roi(5));
maskCS = createMask(roi(6));

roi_cluster = maskV + maskLM + maskAL + maskRL + maskCV + maskCS;
mask_cell = bwlabel(roi_cluster);
figure; imagesc(mask_cell);

ind = find(vasc_mask); mask_cell_V(ind) = 0;
figure; imagesc(mask_cell_V);

%bwout = imCellEditInteractive(max_DF_noVasc);

%% Plot timecourses for events for single days based on ROIs

leverPressTC = stackGetTimeCourses(trialStart_data_avg, mask_cell_AL);
figure; plot(leverPressTC)

sucLeverReleaseTC = stackGetTimeCourses(sucTrialEnd_data_avg, mask_cell_AL);
figure; plot(sucLeverReleaseTC)

earlyLeverReleaseTC = stackGetTimeCourses(earlyTrialEnd_data_avg, mask_cell);
figure; plot(earlyLeverReleaseTC)

missedLeverReleaseTC = stackGetTimeCourses(missedTrialEnd_data_avg, mask_cell);
figure; plot(missedLeverReleaseTC);

sucTrialTargetTC = stackGetTimeCourses(sucTrialTarget_data_avg, mask_cell_AL);
figure; plot(sucTrialTargetTC);

missTrialTargetTC = stackGetTimeCourses(missTrialTarget_data_avg, mask_cell_AL);
figure; plot(missTrialTargetTC);

%% Plot timecourses for events averaging 5 days with mean and sd










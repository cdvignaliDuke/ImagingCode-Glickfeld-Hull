%% get path names
date = '190515';
ImgFolder = strvcat('002');
time = strvcat('1007');
mouse = 'i1304';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

%% load and register
tic
data = [];
CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
%CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
%CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
%CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
load(fName);
nframes = [input.counterValues{end}(end) info.config.frames];
fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
data = sbxread(imgMatFile(1,1:11),0,min(nframes));
if size(data,1)== 2
    data = data(1,:,:,:);
end
data = squeeze(data);
toc

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

data_avg = mean(data(:,:,20001:20500),3);
%% Register data

if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out_bx);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

cStimOne = cell2mat(input.cStimOneOn);
cStimTwo = cell2mat(input.cStimTwoOn);
stimOneContrast = celleqel2mat_padded(input.tStimOneGratingContrast);
stimOneCons = unique(stimOneContrast);
stimOneTime = celleqel2mat_padded(input.tStimOneGratingOnTimeMs);
stimOneTimes = unique(stimOneTime);
nTime = length(stimOneTimes);
nTrials = length(cStimOne);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_one = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & (cStimOne(itrial)+30)<sz(3)
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)-20:cStimOne(itrial)-1),3);
        data_one(:,:,itrial) = mean(data_reg(:,:,cStimOne(itrial)+5:cStimOne(itrial)+25),3);
    end
end
data_one_dfof = (data_one-data_f)./data_f;

if find(stimOneCons ==0)
    zeroConInd = find(stimOneContrast==0);
end

data_dfof_stim = zeros(sz(1),sz(2),nTime+1);
[n n2] = subplotn(nTime+1);
figure;
for it = 1:nTime
    ind = setdiff(find(stimOneTime == stimOneTimes(it)),zeroConInd);
    data_dfof_stim(:,:,it) = nanmean(data_one_dfof(:,:,ind),3);
    subplot(n,n2,it)
    imagesc(data_dfof_stim(:,:,it))
    title(num2str(stimOneTimes(it)))
end
subplot(n,n2,it+1)
ind = setdiff(1:nTrials,zeroConInd);
data_dfof_stim(:,:,it+1) = nanmean(data_one_dfof(:,:,ind),3);
imagesc(data_dfof_stim(:,:,it+1))
title('Avg all')

data_dfof = cat(3,data_dfof_stim(:,:,it+1),data_dfof_stim(:,:,it+1));

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)
% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% Stim two analysis

tc_one = nan(200,nCells,nTrials);
tc_two = nan(30,nCells,nTrials);
   
for itrial = 1:nTrials
    if ~isnan(cStimOne(itrial)) & ~isnan(cStimTwo(itrial)) & (cStimOne(itrial)+189)<sz(3) & (cStimTwo(itrial)+19)<sz(3)
        tc_one(:,:,itrial) = npSub_tc(cStimOne(itrial)-10:cStimOne(itrial)+189,:);
        tc_two(:,:,itrial) = npSub_tc(cStimTwo(itrial)-10:cStimTwo(itrial)+19,:);
    end
end
tc_one_f = mean(tc_one(1:10,:,:));
tc_two_f = mean(tc_two(1:10,:,:));
tc_one_dfof = (tc_one-tc_one_f)./tc_one_f;
tc_two_dfof = (tc_two-tc_two_f)./tc_two_f;

tt1 = [-10:189].*(1000/frame_rate);
tt2 = [-10:19].*(1000/frame_rate);
base_win = [9:14];
resp_win = [16:20];

figure;
tc_one_dfof_all = zeros(200,nCells,nTime);
tc_two_dfof_all = zeros(30,nCells,nTime+1);

for it = 1:nTime
    ind = setdiff(find(stimOneTime == stimOneTimes(it)),zeroConInd);
    tc_one_dfof_all(:,:,it) = nanmean(tc_one_dfof(:,:,ind),3);
    tc_two_dfof_all(:,:,it) = nanmean(tc_two_dfof(:,:,ind),3);
    subplot(2,2,1)
    plot(tt1,squeeze(mean(tc_one_dfof_all(:,:,it)-mean(tc_one_dfof_all(base_win,:,it),1),2)))
    hold on
    subplot(2,2,3)
    plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it)-mean(tc_two_dfof_all(base_win,:,it),1),2)))
    hold on
end
tc_two_dfof_all(:,:,it+1) = nanmean(tc_two_dfof(:,:,zeroConInd),3);
plot(tt2,squeeze(mean(tc_two_dfof_all(:,:,it+1),2)))
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim two')
subplot(2,2,1)
xlabel('Time (ms)')
ylabel('dF/F')
title('Stim one')

tc_one_dfof_resp = mean(tc_one_dfof(resp_win,:,:),1)-mean(tc_one_dfof(base_win,:,:),1);
tc_two_dfof_resp = mean(tc_two_dfof(resp_win,:,:),1)-mean(tc_two_dfof(base_win,:,:),1);
tc_two_dfof_resp_stim = zeros(nCells,nTime+1);
for it = 1:nTime
    ind = setdiff(find(stimOneTime == stimOneTimes(it)),zeroConInd);
    subplot(2,2,2)
    errorbar(stimOneTimes(it),squeeze(mean(nanmean(tc_one_dfof_resp(:,:,ind),3),2)), squeeze(std(nanmean(tc_one_dfof_resp(:,:,ind),3),[],2))./sqrt(nCells),'o')
    hold on
    subplot(2,2,4)
    errorbar(stimOneTimes(it),squeeze(mean(nanmean(tc_two_dfof_resp(:,:,ind),3),2)), squeeze(std(nanmean(tc_two_dfof_resp(:,:,ind),3),[],2))./sqrt(nCells),'o')
    hold on
    tc_two_dfof_resp_stim(:,it) = nanmean(tc_two_dfof_resp(:,:,ind),3);
end
errorbar(0,squeeze(mean(nanmean(tc_two_dfof_resp(:,:,zeroConInd),3),2)), squeeze(std(nanmean(tc_two_dfof_resp(:,:,zeroConInd),3),[],2))./sqrt(nCells),'o')
tc_two_dfof_resp_stim(:,it+1) = nanmean(tc_two_dfof_resp(:,:,zeroConInd),3);
xlabel('Stim one time (ms)')
ylabel('dF/F')
title('Stim two')
xlim([-100 4500])
ylim([0 0.15])
subplot(2,2,2)
xlabel('Stim one time (ms)')
ylabel('dF/F')
title('Stim one')
xlim([-100 4500])
ylim([0 0.15])
suptitle([date ' ' mouse ])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptDuration.pdf']),'-dpdf','-fillpage')

tc_two_dfof_resp_norm = tc_two_dfof_resp_stim./tc_two_dfof_resp_stim(:,end);
figure;
for it = 1:nTime
    errorbar(stimOneTimes(it),squeeze(mean(tc_two_dfof_resp_norm(:,it),1)), squeeze(std(tc_two_dfof_resp_norm(:,it),[],1))./sqrt(nCells),'o')
    hold on
end
xlabel('Stim one time (ms)')
ylabel('Normalized dF/F')
ylim([0 1.2])
set(gca,'xscale','log')

%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'tc_one_dfof', 'tc_two_dfof', 'nCells','frame_rate')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'stimOneTime', 'stimOneTimes', 'nTime', 'stimOneContrast', 'stimOneCons','base_win','resp_win')

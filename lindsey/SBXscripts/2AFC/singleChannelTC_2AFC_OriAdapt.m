%% get path names
date = '181221';
ImgFolder = strvcat('003','004');
time = strvcat('1733');
mouse = 'i1103';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

[s,tUsername] = dos('ECHO %USERNAME%');
if tUsername(1:4) == 'ryan'
    tDir = 'ryan';
elseif tUsername(1:4) == 'lind'
    tDir = 'lindsey';
end
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\' ;
data_base = 'home\lindsey\Data\2P_images\';
%% load dataset
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [fn_base data_base date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ryan\2Pimages\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW14\two-photon imaging\' date '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    if size(time,1) >= irun
        fName = [fn_base 'Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
        load(fName);
        temp(irun) = input;
        if irun>1
            if isfield(input, 'tLeftTrial')
                ntrials = size(input.trialOutcomeCell,2);
                for itrial = 1:ntrials
                    temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                    temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset;
                    temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                    temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset;
                end
            end
        end
    end
    
    nframes = min([info.config.frames input.counterValues{end}(end)]);
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    toc
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Plot outcome by trial number
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');

figure;
plot(smooth(SIx,10));
hold on
plot(smooth(MIx,10));

%% Crop data and input struct
%change trial range in trialChopper
input = trialChopper(input,[1 415]);
data = data(:,:,input.counterValues{1}(1):input.counterValues{end}(end));

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%% Register data

data_avg = mean(data(:,:,20001:20500),3);

[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data

%% test stability

figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
%% find activated cells

tGratingOri = celleqel2mat_padded(input.tGratingDirectionStart);
Oris = unique(tGratingOri);
nOri = length(Oris);
cStimOn = celleqel2mat_padded(input.cStimOn);
nTrials = length(tGratingOri);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3);
    data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3);
end
data_targ_dfof = (data_targ-data_f)./data_f;
[n n2] = subplotn(nOri+1);
figure;
data_dfof = zeros(sz(1),sz(2),nOri+1);
for iori = 1:nOri
    subplot(n,n2,iori)
    ind = find(tGratingOri == Oris(iori));
    data_dfof(:,:,iori)= mean(data_targ_dfof(:,:,ind),3);
    imagesc(data_dfof(:,:,iori));
    title([num2str(Oris(iori)) ' deg'])
end
subplot(n,n2,iori+1)
data_dfof(:,:,iori+1)= mean(data_targ_dfof,3);
imagesc(data_dfof(:,:,iori+1));
title('All')

data_dfof_max = max(data_dfof,[],3);
figure; imagesc(data_dfof_max)
data_dfof = cat3(data_dfof_max,data_dfof);
    
%% cell segmentation 

mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell)

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
%np_tc = stackGetTimeCourses(data_reg,mask_np);
clear np_tc np_tc_down
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
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

save(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_tc np_tc

%% 2AFC analysis
frameRateHz = double(input.frameRateHz);
cStimOn = celleqel2mat_padded(input.cStimOn);
tGratingOri = celleqel2mat_padded(input.tGratingDirectionStart);
tOris = unique(tGratingOri);
nOri = length(tOris);
aGratingOri = celleqel2mat_padded(input.aGratingDirectionDeg);
aGratingContrast = celleqel2mat_padded(input.aGratingContrast);
aCons = unique(aGratingContrast);
naCon = length(aCons);
aOris = unique(aGratingOri);
naOri = length(aOris);
nCells = size(npSub_tc,2);
nframes = size(npSub_tc,1);
nTrials = size(aGratingOri,2);
data_stim = nan(50,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
    end
end
dataf = mean(data_stim(1:20,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
tt = [-20:29].*(1000./frameRateHz);
figure;
plot(nanmean(mean(data_stim_dfof,2),3));
vline([16 20 25 29])

base_win = [16:20];
resp_win = [25:29];

nt = cell(3,nOri);
x = 1;
start = 1;
figure;
for icon = 1:naCon
    if aCons(icon) == 1
        for iaOri = 1:naOri
            ind_con = find(aGratingContrast == aCons(icon));
            ind_con = intersect(ind_con, find(aGratingOri == aOris(iaOri)));
            for iori = 1:nOri
                ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
                subplot(3,nOri,start)
                plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
                hold on
                title(['Adapt: Ori = '  num2str(aOris(iaOri))])
                ylim([-0.01 0.1])
                start = start+1;
                nt{x,iori} = ind_ori;
            end
            x = 1+x;
        end
    else
        ind_con = find(aGratingContrast == aCons(icon));
        for iori = 1:nOri
            ind_ori = intersect(ind_con, find(tGratingOri == tOris(iori)));
            subplot(3,nOri,start)
            plot(tt, mean(nanmean(data_stim_dfof(:,:,ind_ori),3),2))
            hold on
            title(['No adapt; Ori = ' num2str(tOris(iori))])
            ylim([-0.01 0.1])
            start = start+1;
            nt{x,iori} = ind_ori;
        end
        x = 1+x;
    end
    legend({'No adapt','0deg','90deg'})
end
print(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse_byTarget.pdf']), '-dpdf','-bestfit')

ind_aCon0 = find(aGratingContrast == 0);
ind_cond{1} = ind_aCon0;
ind_cond{2} = intersect(find(aGratingContrast),find(aGratingOri==0));
ind_cond{3} = intersect(find(aGratingContrast),find(aGratingOri==90));

h_ori = nan(nOri, nCells,3);
[h_all p_all] = ttest(mean(data_stim_dfof(base_win,:,ind_cond{i}),1),mean(data_stim_dfof(resp_win,:,ind_cond{i}),1),'tail','left','dim',3);
for iori = 1:nOri
    ind_ori = intersect(ind_cond{1}, find(tGratingOri == tOris(iori))); 
    [h_ori(iori,:) p_ori(iori,:)] = ttest(mean(data_stim_dfof(base_win,:,ind_ori),1),mean(data_stim_dfof(resp_win,:,ind_ori),1),'tail','left','dim',3,'alpha',0.05./(nOri-1));
end

good_ind = find(sum([h_all; h_ori],1));
data_stim_resp = squeeze(mean(data_stim_dfof(resp_win,:,:),1)-mean(data_stim_dfof(base_win,:,:),1));

data_resp_ori = zeros(3,nOri,length(good_ind));

[n n2] = subplotn(length(good_ind));
h_anova = nan(1,nCells);
figure;
for iCell = 1:length(good_ind)
    subplot(n,n2,iCell)
    iC = good_ind(iCell);
    temp_resp = [];
    for i = 1:3
        for iOri = 1:nOri
            data_resp_ori(i,iOri,iCell) = mean(data_stim_resp(iC,nt{i,iOri}),2);
            if i == 1
                temp_resp = [temp_resp; data_stim_resp(iC,nt{i,iOri})' iOri.*ones(length(nt{i,iOri}),1)];
            end
        end
        h_anova(iCell) = anova1(temp_resp(:,1), temp_resp(:,2),'off');
        plot(tOris,data_resp_ori(i,:,iCell),'-o')
        hold on
    end
    ylim([-0.1 1])
end

[ori_max ori_pref] = max(data_resp_ori(1,:,:),[],2);
figure
for ipref = 1:nOri
    ind = intersect(find(h_anova(good_ind)),find(squeeze(ori_pref) == ipref));
    subplot(3,3,ipref)
    for i = 1:3
        plot(tOris, mean(data_resp_ori(i,:,ind),3),'-o')
        hold on
    end
    title([num2str(tOris(ipref)) 'deg- n = ' num2str(length(ind)) ' cells'])
end

aCon_p1 = [NaN aGratingContrast];
aOri_p1 = [NaN aGratingOri];
ind_aCon0_p1{1} = intersect(ind_aCon0,find(aCon_p1==0));
ind_aCon0_p1{2} = intersect(ind_aCon0,intersect(find(aCon_p1==1),find(aOri_p1 ==0)));
ind_aCon0_p1{3} = intersect(ind_aCon0,intersect(find(aCon_p1==1),find(aOri_p1 ==90)));

data_aCon0_ori = zeros(3,nOri,length(good_ind));
figure;
for iCell = 1:length(good_ind)
    subplot(n,n2,iCell)
    iC = good_ind(iCell);
    for i = 1:3
        for iOri = 1:nOri
            data_aCon0_ori(i,iOri,iCell) = mean(data_stim_resp(iC,intersect(ind_aCon0_p1{i},nt{1,iOri})),2);
        end
        plot(tOris,data_aCon0_ori(i,:,iCell),'-o')
        hold on
    end
end

figure;
for iCell = 1:length(good_ind)
    subplot(n,n2,iCell)
    iC = good_ind(iCell);
    for i = 1:3
        plot(tt, mean(data_stim_dfof(:,iC,ind_aCon0_p1{i}),3))
        hold on
    end
end

figure;
subplot(1,2,1)
for i = 1:3
 plot(tt, mean(mean(data_stim_dfof(:,:,ind_cond{i}),3),2))
 hold on
end
title('Current trial')
ylabel('dF/F')
xlabel('Time from target (ms)')
legend({'Con = 0','Ori = 0', 'Ori = 90'})
subplot(1,2,2)
for i = 1:3
 plot(tt, mean(mean(data_stim_dfof(:,:,ind_aCon0_p1{i}),3),2))
 hold on
end
title('Previous trial, for Current Con = 0')
ylabel('dF/F')
xlabel('Time from target (ms)')
print(fullfile([fn_base 'home\lindsey\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_avgTargetResponse.pdf']), '-dpdf','-bestfit')



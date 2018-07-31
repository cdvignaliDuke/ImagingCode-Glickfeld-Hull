%expt info
date = '170704';
mouse = 'i720';
run_str = 'runs-003-004';
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% extract stim parameters
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
if iscell(input.nFramesOn)
    nOn = input.nFramesOn{1};
else
    nOn = input.nFramesOn;
end

tCyc = cell2mat(input.tCyclesOn);
cStart = celleqel2mat_padded(input.cFirstStim);
cTarget = celleqel2mat_padded(input.cTargetOn);
nTrials = length(tCyc);
nCells = size(npSub_tc,2);
maxCyc = max(tCyc,[],2);
tFramesOff = nan(nTrials,maxCyc);
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');
FIx = strcmp(input.trialOutcomeCell, 'failure');
nCyc = tCyc;
nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
for itrial = 1:nTrials
    if isfield(input, 'tFramesOff')
        if length(input.tFramesOff{itrial}>0)
            tempFramesOff = input.tFramesOff{itrial};
        else
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
            input.tFramesOff{itrial} = tempFramesOff;
        end
    else
        if iscell(input.nFramesOff)
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
        else
            tempFramesOff = input.nFramesOff.*(ones(1,tCyc(itrial)));
        end
    end

    tFramesOff(itrial,1:tCyc(itrial)) = tempFramesOff(1:tCyc(itrial));
end
targCon = celleqel2mat_padded(input.tGratingContrast);
if isfield(input,'doRandCon') & input.doRandCon
    baseCon = nan(maxCyc,nTrials);
    for itrial = 1:nTrials
        baseCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
    end
    ind_con = [];
else
    baseCon = celleqel2mat_padded(input.tBaseGratingContrast);
    ind_con = intersect(find(targCon == 1),find(baseCon == 0));
end
baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
tGratingDir = round(double(celleqel2mat_padded(input.tGratingDirectionDeg)),0);
if sum(tGratingDir-baseDir) == 0
    targetDelta = tGratingDir-baseDir;
else
    targetDelta = tGratingDir;
end
deltas = unique(targetDelta);
nDelta = length(deltas);
offs = unique(tFramesOff(:,1));
noff = length(offs);
frameRateHz = input.frameRateHz;

%create trial dfof from params
data_trial = nan(120,nCells,maxCyc+1,nTrials);
for itrial = 1:nTrials
    tempFramesOff = tFramesOff(itrial,1:tCyc(itrial));
    if ~isnan(cStart(itrial))
        for icyc = 1:nCyc(itrial)
            if icyc > 1
                cyc_add = ((icyc-1)*nOn)+sum(tempFramesOff(1:icyc-1));
            else
                cyc_add = 0;
            end
            if cStart(itrial)+99+cyc_add <= size(npSub_tc,1)
                data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-20+cyc_add:cStart(itrial)+99+cyc_add,:);
            else
                data_trial(:,:,icyc,itrial) = NaN(120,nCells);
            end 
        end
    else
        data_trial(:,:,icyc,itrial) = NaN(120,nCells);
    end
end
data_f = nanmean(data_trial(1:20,:,1,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

%these need to be adjusted for the experiment
base_win =21:23;
resp_win =27:29;

figure;
subplot(2,1,1)
plot(squeeze(nanmean(mean(data_dfof(:,:,1,:),2),4)));
vline(base_win)
vline(resp_win)
title('Baseline')
subplot(2,1,2)
sz = size(data_dfof);
data_targ = zeros(sz(1),sz(2),length([find(SIx)]));
for itrial = 1:sz(4);
    if find([find(SIx)] == itrial)
        data_targ(:,:,itrial) = data_dfof(:,:,nCyc(itrial),itrial);
    end
end
plot(squeeze(nanmean(mean(data_targ,2),3)));
title('Target')
vline(base_win)
vline(resp_win)

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'baseCon', 'ind_con', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart', 'base_win','resp_win')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof')

    clear data_tc data_f data_targ np_tc npSub_tc data_trial
%% Find responsive cells
%base responsive cells
[x1,y1] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
good_ind_temp = find(x1);

%find late responding cells and remove
tc_all = squeeze(nanmean(mean(bsxfun(@minus,data_dfof(:,:,1,:),mean(data_dfof(base_win,:,1,:),1)),3),4));
resp_diff = diff(tc_all);
[max_val max_time] = max(resp_diff(20:end,:),[],1);

ind1 = find(max_time<base_win(end)-20);
ind2 = find(max_time>resp_win(end)-20);
ind3 = setdiff(1:nCells, [ind1 ind2]);
good_ind_base = intersect(good_ind_temp,ind3);

figure;
subplot(2,2,1)
plot(tc_all(:,ind1))
subplot(2,2,2)
plot(tc_all(:,ind2))
subplot(2,2,3)
plot(tc_all(:,ind3))
subplot(2,2,4)
plot(tc_all(:,good_ind_base))

%find target responsive cells
[x2,y2] = ttest(squeeze(mean(data_dfof(base_win,:,6,:),1))',squeeze(mean(data_dfof(resp_win,:,6,:),1))','tail','left','alpha',0.05);

h1 = zeros(ndir,nCells);
p1 = zeros(ndir,nCells);
data_dfof_delta = nan(size(data_dfof,1), nCells,nDelta);
for idelta = 1:nDelta
    ind = find(targetDelta == deltas(idelta));
    data_dfof_delta(:,:,idelta) = nanmean(data_dfof(:,:,2,ind),4);
    for iCell = 1:nCells
        [h1(idelta,iCell), p1(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,2,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,2,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
    end
end
good_ind_temp = find(x2+sum(h1,1));

%find late responding cells and remove
tc_all = squeeze(nanmean(mean(bsxfun(@minus,data_dfof(:,:,6,:),mean(data_dfof(base_win,:,6,:),1)),3),4));
resp_diff = diff(tc_all);
[max_val max_time] = max(resp_diff(20:end,:),[],1);

ind1 = find(max_time<base_win(end)-20);
ind2 = find(max_time>resp_win(end)-20);
ind3 = setdiff(1:nCells, [ind1 ind2]);
good_ind_targ = intersect(good_ind_temp,ind3);

figure;
subplot(2,2,1)
plot(tc_all(:,ind1))
subplot(2,2,2)
plot(tc_all(:,ind2))
subplot(2,2,3)
plot(tc_all(:,ind3))
subplot(2,2,4)
plot(tc_all(:,good_ind_targ))

%% auROC of baseline stimuli for n interval

p1_resp = squeeze(nanmean(mean(data_dfof(resp_win,:,1,:),1)-mean(data_dfof(base_win,:,1,:),1),4));

for ioff = 1:noff
    resp{ioff} = [];
    base{ioff} = [];
    for icyc = 2:maxCyc-1
        ind = find(tFramesOff(:,icyc) == offs(ioff));
        resp{ioff} = cat(3, resp{ioff}, squeeze(data_dfof(:,:,icyc+1,ind)));
        base{ioff} = cat(3, base{ioff}, squeeze(data_dfof(:,:,icyc,ind)));
    end
end

roc_base = nan(noff,nCells);
for ioff = 1:noff
    resp_temp = resp{ioff};
    base_temp = base{ioff};
    resp_diff = squeeze(mean(resp_temp(resp_win,:,:),1)-mean(resp_temp(base_win,:,:),1));
    base_diff = squeeze(mean(base_temp(resp_win,:,:),1)-mean(base_temp(base_win,:,:),1));
    for iCell = 1:length(good_ind_base)
        iC = good_ind_base(iCell);
        roc_base(ioff,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
    end
end

figure; 
errorbar(offs*(1000/frameRateHz), nanmean(roc_base,2), nanstd(roc_base,[],2)./sqrt(length(good_ind_base)),'-ok');
title([mouse ' ' date ' ' run_str ' - Base responses- All n-1'])
xlabel('N interval')
ylim([0.4 0.6])

print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Base_allN-1.pdf']),'-dpdf','-fillpage')

%% auROC of baseline stimuli for all n interval
for ioff = 1:noff
    resp{ioff} = [];
    base{ioff} = [];
    for icyc = 2:maxCyc-1
        ind = find(tFramesOff(:,icyc-1) == offs(ioff));
        resp{ioff} = cat(3, resp{ioff}, squeeze(data_dfof(:,:,icyc+1,ind)));
        base{ioff} = cat(3, base{ioff}, squeeze(data_dfof(:,:,icyc,ind)));
    end
end

roc_base_allN = nan(noff,nCells);
for ioff = 1:noff
    resp_temp = resp{ioff};
    base_temp = base{ioff};
    resp_diff = squeeze(mean(resp_temp(resp_win,:,:),1)-mean(resp_temp(base_win,:,:),1));
    base_diff = squeeze(mean(base_temp(resp_win,:,:),1)-mean(base_temp(base_win,:,:),1));
    for iCell = 1:length(good_ind_base)
        iC = good_ind_base(iCell);
        roc_base_allN(ioff,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
    end
end

figure; 
errorbar(offs*(1000/frameRateHz), nanmean(roc_base_allN,2), nanstd(roc_base_allN,[],2)./sqrt(length(good_ind_base)),'-ok');
title([mouse ' ' date ' ' run_str ' - Base responses- All N'])
xlabel('N-1 interval')
ylim([0.4 0.6])

print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Base_allN.pdf']),'-dpdf','-fillpage')
%% auROC of baseline stimuli for n-1 interval

resp = cell(noff, noff);
base = cell(noff, noff);
for ioff_N = 1:noff
    for ioff_N1 = 1:noff
        for icyc = 2:maxCyc-1
            ind = intersect(find(tFramesOff(:,icyc) == offs(ioff_N)),find(tFramesOff(:,icyc-1) == offs(ioff_N1)));
            resp{ioff_N, ioff_N1} = cat(3, resp{ioff_N, ioff_N1}, squeeze(data_dfof(:,:,icyc+1,ind)));
            base{ioff_N, ioff_N1} = cat(3, base{ioff_N, ioff_N1}, squeeze(data_dfof(:,:,icyc,ind)));
        end
    end
end

roc_base_N1 = nan(noff,noff,nCells);
base_resp = nan(noff,noff,nCells);
base_resp_N1 = nan(noff,noff,nCells);
fs_baseResp = cell(noff,noff);
for ioff_N = 1:noff
    for ioff_N1 = 1:noff
        resp_temp = resp{ioff_N, ioff_N1};
        base_temp = base{ioff_N, ioff_N1};
        resp_diff = squeeze(mean(resp_temp(resp_win,:,:),1)-mean(resp_temp(base_win,:,:),1));
        base_diff = squeeze(mean(base_temp(resp_win,:,:),1)-mean(base_temp(base_win,:,:),1));
        for iCell = 1:length(good_ind_base)
            iC = good_ind_base(iCell);
            roc_base_N1(ioff_N,ioff_N1,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
        end
        base_resp(ioff_N,ioff_N1,:) = reshape(mean(resp_diff,2)./p1_resp', [1 1 nCells]);
        base_resp_N1(ioff_N,ioff_N1,:) = reshape(mean(base_diff,2)./p1_resp', [1 1 nCells]);
        fs_baseResp{ioff_N,ioff_N1} = cat(3,base_diff,resp_diff);
    end
end

figure; 
for ioff_N = 1:noff
    subplot(1,noff,ioff_N)
    errorbar(offs*(1000/frameRateHz), nanmean(roc_base_N1(ioff_N,:,:),3), nanstd(roc_base_N1(ioff_N,:,:),[],3)./sqrt(length(good_ind_base)),'-ok');
    title(['N = ' num2str(offs(ioff_N)*(1000/frameRateHz))])
    xlabel('N-1 interval')
    ylim([0.4 0.6])
end
suptitle([mouse ' ' date ' ' run_str '- Base responses'])
print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Base_byN-1.pdf']),'-dpdf','-fillpage')

%% auROC of target stimuli for n interval

roc_targ = nan(noff,nDelta,nCells);
for ioff = 1:noff
    for idelta = 1:nDelta
        ind = intersect(find(tGratingDir == deltas(idelta)), find(tFramesOff(:,5) == offs(ioff)));
        resp_diff = squeeze(mean(data_dfof(resp_win,:,6,ind),1)- mean(data_dfof(base_win,:,6,ind),1));
        base_diff = squeeze(mean(data_dfof(resp_win,:,5,ind),1)- mean(data_dfof(base_win,:,5,ind),1));
        for iCell = 1:length(good_ind_targ)
            iC = good_ind_targ(iCell);
            roc_targ(ioff,idelta,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
        end
    end
end

figure;
for idelta = 1:nDelta
    subplot(1,2,idelta)
    errorbar(offs*(1000/frameRateHz), nanmean(roc_targ(:,idelta,:),3), nanstd(roc_targ(:,idelta,:),[],3)./sqrt(length(good_ind_targ)),'-ok');
    title([num2str(deltas(idelta)) ' deg'])
    xlabel('N interval')
    ylim([0.3 0.7])
end
suptitle([mouse ' ' date ' ' run_str ' - Target responses- All n-1'])
print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Targ_allN-1.pdf']),'-dpdf','-fillpage')

%% auROC of target stimuli for all n interval

roc_targ_allN = nan(noff,nDelta,nCells);
for ioff = 1:noff
    for idelta = 1:nDelta
        ind = intersect(find(tGratingDir == deltas(idelta)), find(tFramesOff(:,4) == offs(ioff)));
        resp_diff = squeeze(mean(data_dfof(resp_win,:,6,ind),1)- mean(data_dfof(base_win,:,6,ind),1));
        base_diff = squeeze(mean(data_dfof(resp_win,:,5,ind),1)- mean(data_dfof(base_win,:,5,ind),1));
        for iCell = 1:length(good_ind_targ)
            iC = good_ind_targ(iCell);
            roc_targ_allN(ioff,idelta,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
        end
    end
end

figure;
for idelta = 1:nDelta
    subplot(1,2,idelta)
    errorbar(offs*(1000/frameRateHz), nanmean(roc_targ_allN(:,idelta,:),3), nanstd(roc_targ_allN(:,idelta,:),[],3)./sqrt(length(good_ind_targ)),'-ok');
    title([num2str(deltas(idelta)) ' deg'])
    xlabel('N interval')
    ylim([0.3 0.7])
end
suptitle([mouse ' ' date ' ' run_str ' - Target responses- All n'])
print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Targ_allN.pdf']),'-dpdf','-fillpage')

%% auROC of target stimuli for n-1 interval

roc_targ_N1 = nan(noff,noff,nDelta,nCells);
targ_resp = nan(noff,noff,nDelta,nCells);
targ_resp_N1 = nan(noff,noff,nDelta,nCells);
fs_targResp = cell(noff,noff,nDelta);
for ioff_N = 1:noff
    for ioff_N1 = 1:noff
        for idelta = 1:nDelta
            ind = intersect(find(tGratingDir == deltas(idelta)), intersect(find(tFramesOff(:,5) == offs(ioff_N)), find(tFramesOff(:,4) == offs(ioff_N1))));
            resp_diff = squeeze(mean(data_dfof(resp_win,:,6,ind),1)- mean(data_dfof(base_win,:,6,ind),1));
            base_diff = squeeze(mean(data_dfof(resp_win,:,5,ind),1)- mean(data_dfof(base_win,:,5,ind),1));
            for iCell = 1:length(good_ind_targ)
                iC = good_ind_targ(iCell);
                roc_targ_N1(ioff_N,ioff_N1,idelta,iC) = roc_gh(base_diff(iC,:), resp_diff(iC,:));
            end
            targ_resp(ioff_N,ioff_N1,idelta,:) = reshape(mean(resp_diff,2)./p1_resp',[1 1 nCells]);
            targ_resp_N1(ioff_N,ioff_N1,idelta,:) = reshape(mean(base_diff,2)./p1_resp',[1 1 nCells]);
            fs_targResp{ioff_N,ioff_N1,idelta} = cat(3, base_diff, resp_diff); 
        end
    end
end

figure;
i = 1;
for idelta = 1:nDelta
    for ioff_N = 1:ioff
        subplot(2,3,i)
        errorbar(offs*(1000/frameRateHz), squeeze(nanmean(roc_targ_N1(ioff_N,:,idelta,:),4)), squeeze(nanstd(roc_targ_N1(ioff_N,:,idelta,:),[],4)./sqrt(length(good_ind_targ))),'-ok');
        title([num2str(deltas(idelta)) ' deg'])
        title([num2str(deltas(idelta)) ' deg- N = ' num2str(offs(ioff_N)*(1000/frameRateHz))])
        xlabel('N-1 interval')
        ylim([0.3 0.7])
        i = 1+i;
    end
end

suptitle([mouse ' ' date ' ' run_str '- Target responses'])

print(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Targ_byN-1.pdf']),'-dpdf','-fillpage')

save(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ROC.mat']),'base_resp','targ_resp','base_resp_N1','targ_resp_N1','roc_base','roc_base_N1','roc_base_allN','roc_targ','roc_targ_N1','roc_targ_allN','good_ind_targ', 'good_ind_base')


%% average response by cycle number by interval
resp_dfof = squeeze(mean(data_dfof(resp_win,:,:,:),1)-mean(data_dfof(base_win,:,:,:),1));

tFramesOff_ioff = cell(4,noff);
ntrials = zeros(4,noff);
resp_dfof_ioff = nan(nCells,5,noff);
for icyc = 1:4
    for ioff = 1:noff
        tFramesOff_ioff{icyc,ioff} = find(sum(tFramesOff(:,1:icyc) == offs(ioff),2) == icyc);
        ntrials(icyc,ioff) = length(tFramesOff_ioff{icyc,ioff});
        resp_dfof_ioff(:,icyc+1,ioff) = nanmean(resp_dfof(:,icyc+1,tFramesOff_ioff{icyc,ioff}),3);
        if icyc == 1
            resp_dfof_ioff(:,1,ioff) = nanmean(resp_dfof(:,1,:),3);
        end
    end
end

baseresp_dfof_norm = bsxfun(@rdivide, resp_dfof(:,1:5,:), mean(resp_dfof(:,1,:),3));
figure; 
errorbar(1:5, mean(nanmean(baseresp_dfof_norm(good_ind_base,:,:),3),1), std(nanmean(baseresp_dfof_norm(good_ind_base,:,:),3),[],1)./sqrt(length(good_ind_base)),'-o')
ylim([0 1])

baseresp_dfof_ioff_norm = bsxfun(@rdivide, resp_dfof_ioff, nanmean(resp_dfof_ioff(:,1,:),3));
figure; 
for ioff = 1:noff
    errorbar(1:5, mean(baseresp_dfof_ioff_norm(good_ind_base, :,ioff),1), std(baseresp_dfof_ioff_norm(good_ind_base, :,ioff),[],1)./sqrt(length(good_ind_base)),'-o')
    hold on
end
ylim([0 1])

temp_off = zeros(size(tFramesOff(:,1:4)));
for ioff = 1:noff
    ind = find(tFramesOff(:,1:4) == offs(ioff));
    temp_off(ind) = ioff;
end
tTrialLength = (temp_off.*250)+100;
tTotalTrialLength = cumsum(tTrialLength,2);

edges = [0 350:500:4000];
sz = length(edges);
[n bin] = histc(tTotalTrialLength,edges);
baseresp_dfof_ibin_norm = nan(nCells,length(n)-1);
for ibin = 1:length(n)-1
    ind = [];
    if ibin ==1
        baseresp_dfof_ibin_norm(:,1) = mean(baseresp_dfof_norm(:,1,:),3);
    else
        for icyc = 1:4
            ind = [ind; find(bin(:,icyc) == ibin)];
        end
        baseresp_dfof_ibin_norm(:,ibin) = mean(baseresp_dfof_norm(:,icyc,ind),3);
    end
end

figure; 
x = [0 mean([edges(2:sz-1);edges(3:end)],1)];
errorbar(x, mean(baseresp_dfof_ibin_norm(good_ind_base,:),1), std(baseresp_dfof_ibin_norm(good_ind_base,:),[],1)./sqrt(length(good_ind_base)),'-o')
ylim([0 1])

baseresp_dfof_norm = nanmean(bsxfun(@rdivide, resp_dfof(:,1:5,:), mean(resp_dfof(:,1,:),3)),3);
save(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_baseResp.mat']),'baseresp_dfof_norm','baseresp_dfof_ioff_norm','baseresp_dfof_ibin_norm','x', 'edges')

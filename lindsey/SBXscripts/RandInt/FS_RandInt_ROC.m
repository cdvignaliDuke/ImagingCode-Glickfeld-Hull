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

print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Base_allN-1.pdf']),'-dpdf','-fillpage')

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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Base_byN-1.pdf']),'-dpdf','-fillpage')

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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Targ_allN-1.pdf']),'-dpdf','-fillpage')


%% auROC of target stimuli for n-1 interval

roc_targ_N1 = nan(noff,noff,nDelta,nCells);
targ_resp = nan(noff,noff,nDelta,nCells);
targ_resp_N1 = nan(noff,noff,nDelta,nCells);
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

print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auROC_Targ_byN-1.pdf']),'-dpdf','-fillpage')

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ROC.mat']),'base_resp','targ_resp','base_resp_N1','targ_resp_N1','roc_base','roc_base_N1','roc_targ','roc_targ_N1', 'good_ind_targ', 'good_ind_base')



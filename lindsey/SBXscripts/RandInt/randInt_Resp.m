%find good cells- responsive to all base or at least one base direction
[x,y] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
h = zeros(ndir,nCells);
p = zeros(ndir,nCells);


if ndir > 1
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
        data_dfof_dir = nan(40, nCells, maxCyc, ndir);
        for icyc = 1:maxCyc
            ind2 = intersect(ind,find(nCyc >= icyc));
            data_dfof_dir(:,:,icyc,idir)= nanmean(nanmean(data_dfof(:,:,icyc,ind2),3),4);
        end
        for iCell = 1:nCells
            [h(idir,iCell), p(idir,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,1,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,1,ind),1)),'tail','left','alpha',0.05./(ndir-1));
        end
    end
    good_ind = unique([find(x) find(sum(h,1)>0)]);
else
    data_dfof_dir = nan(40, nCells, maxCyc);
    ind_cyc = zeros(1,maxCyc);
    for icyc = 1:maxCyc
        if icyc == 1
            %ind = intersect(find(SIx),find(nCyc>icyc));
            ind = find(tCyc>icyc);
            ind_cyc(1,icyc)= length(ind);
        else
            %ind = intersect(find(SIx), find(nCyc >= icyc));
            ind = find(tCyc >= icyc);
            ind_cyc(1,icyc)= length(ind);
        end
        data_dfof_dir(:,:,icyc)= squeeze(nanmean(data_dfof(:,:,icyc,ind),4));
    end
    good_ind = find(x);
end

%plot ori tuning of all good cells
if ndir > 1
figure; 
[n n2] = subplotn(length(good_ind));
data_dir_avg = zeros(nCells,ndir);
max_dir = NaN(nCells,1);
for i = 1:length(good_ind)
    subplot(n,n2,i); 
    iC = good_ind(i);
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
        data_dir_avg(iC,idir) = squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1));
        if h(idir,iC)
            errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'or');
        else
            errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'ok');
        end
        hold on
        h_ind = find(h(:,iC));
        if length(h_ind)>0
            [max_val max_ind] = max(data_dir_avg(iC,h_ind),[],2);
            max_dir(iC,:) = h_ind(max_ind);
        else
            [max_val max_ind] = max(data_dir_avg(iC,:),[],2);
            max_dir(iC,:) = max_ind;
        end
    end
    title(['Cell ' num2str(iC) '- Pref ' num2str(dirs(max_dir(iC,:)))])
    ylim([0 max(mean(data_dfof_dir(31:35,iC,:),1),[],3)*1.5])
    xlim([-dirs(2) dirs(end)+dirs(2)])
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
elseif ~doFromRef
max_dir = ones(nCells,1);
end
%response by cycle- avg all cells
tt = [-20:19].*frameRateHz;
data_dfof_base = zeros(40,nCells,maxCyc);
if ~doFromRef
for i = 1:length(good_ind)
    iC = good_ind(i);
    data_dfof_base(:,iC,:) = (data_dfof_dir(:,iC,:,max_dir(iC,:)));
end
else
%data_dfof_base = data_dfof_dir(:,:,1:maxCyc+1);
end

figure;
subplot(2,1,1)
for icyc = 1:maxCyc
plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,icyc),nanmean(data_dfof_base(base_win,good_ind,icyc),1)),2))
hold on
end
ylim([-0.02 0.15])
xlabel('Time (ms)')
legend([repmat('cycle ',[maxCyc 1]) num2str([1:maxCyc]')])
legend('Location','Northwest')
subplot(2,1,2)
data_dfof_base_diff = zeros(maxCyc,length(good_ind));
for icyc = 1:maxCyc
data_dfof_base_diff(icyc,:) = nanmean(data_dfof_base(resp_win,good_ind,icyc),1)- nanmean(data_dfof_base(base_win,good_ind,icyc),1);
errorbar(icyc,nanmean(data_dfof_base_diff(icyc,:),2),nanstd(data_dfof_base_diff(icyc,:),[],2)./sqrt(length(good_ind)),'o')
hold on
end
ylim([0 0.15])
ylabel('dF/F')
xlabel('cycle #')
set(gca,'XTick',1:maxCyc)
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells- ' num2str(ind_cyc)]) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc.pdf']),'-dpdf','-bestfit')

figure
[n n2] = subplotn(length(good_ind));
for iCell = 1:length(good_ind)
iC = good_ind(iCell);
subplot(n,n2,iCell)
for icyc = 1:2
    plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,iC,icyc),nanmean(data_dfof_base(base_win,iC,icyc),1)),2))
    hold on
end
end
suptitle([mouse ' ' date])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_2Cyc_byCell.pdf']),'-dpdf','-bestfit')

off_str{1} = 'base';
for ioff = 1:noff
off_str{ioff+1} = num2str(offs(ioff).*frameRateHz);
end
%response by off frames, avg all cycles all cells
data_dfof_off = zeros(40,nCells,noff);
for i = 1:length(good_ind)
iC = good_ind(i);
for ioff = 1:noff
    temp = [];
    for icyc = 1:maxCyc-1
        if ~doFromRef
            ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
        else
            ind = find(tFramesOff(:,icyc) == offs(ioff));
        end
        temp = [temp squeeze(data_dfof(:,iC,icyc+1,ind))];
    end
    data_dfof_off(:,iC,ioff) = nanmean(temp,2);
end
end
figure;
subplot(2,1,1)
plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
hold on
for ioff = 1:noff
plot(tt,nanmean(bsxfun(@minus,data_dfof_off(:,good_ind,ioff), nanmean(data_dfof_off(base_win,good_ind,ioff),1)),2));
end
legend(off_str)
legend('Location','Northwest')
ylim([-0.02 0.15])
xlabel('Time (ms)')
subplot(2,1,2)
errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
hold on
data_dfof_base_off = zeros(noff,length(good_ind));
for ioff = 1:noff
temp = nanmean(data_dfof_off(28:31,good_ind,ioff),1)- nanmean(data_dfof_off(20:23,good_ind,ioff),1);
errorbar((offs(ioff)*frameRateHz)./1000,mean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
hold on
end
set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
xlim([0.05 10])
ylim([0 0.15])
ylabel('dF/F')
xlabel('Off Interval (ms)')
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt.pdf']),'-dpdf', '-bestfit')

%response by off frames by cell, avg all cycles
figure;
[n n2] = subplotn(length(good_ind));
for iCell = 1:length(good_ind)
iC = good_ind(iCell);
subplot(n, n2, iCell)
plot(tt,bsxfun(@minus,data_dfof_base(:,iC,1),nanmean(data_dfof_base(base_win,iC,1),1)))
hold on
for ioff = 1:noff
    plot(tt,bsxfun(@minus,data_dfof_off(:,iC,ioff), nanmean(data_dfof_off(base_win,iC,ioff),1)));
end
end
suptitle([mouse ' ' date]) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_byCell.pdf']),'-dpdf', '-bestfit')


%response by off frames by cycle, avg all cells
data_dfof_off_cyc = zeros(40,nCells,noff,maxCyc);
ind_n = zeros(maxCyc,noff);
for i = 1:length(good_ind)
iC = good_ind(i);
for ioff = 1:noff
    for icyc = 1:maxCyc
        if ~doFromRef
            ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
            ind_n(icyc,ioff) = length(ind);
        else
            ind = find(tFramesOff(:,icyc) == offs(ioff));
        end
        data_dfof_off_cyc(:,iC,ioff,icyc) = nanmean(data_dfof(:,iC,icyc+1,ind),4);
    end
end
end
figure;
[n n2] = subplotn(double(maxCyc-1));
for icyc = 1:maxCyc-1
subplot(n,n2,double(icyc))
plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
hold on
for ioff = 1:noff
    plot(tt,mean(bsxfun(@minus,data_dfof_off_cyc(:,good_ind,ioff,icyc),mean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1)),2));
end
if icyc<(maxCyc-1)
legend(off_str,'Location','Northwest')
end
title(['Cycle ' num2str(icyc+1) '- ' num2str(ind_n(icyc,:)) ' trials'])
ylim([-0.02 0.15])
xlabel('Time (ms)')
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt.pdf']),'-dpdf', '-bestfit')

figure;
for icyc = 1:(maxCyc-1)
subplot(n,n2,double(icyc))
errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
hold on
for ioff = 1:noff
    temp = bsxfun(@minus, nanmean(data_dfof_off_cyc(resp_win,good_ind,ioff,icyc),1), nanmean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1));
    errorbar((offs(ioff)*frameRateHz)./1000,nanmean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
    hold on
end
set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
xlim([0.05 10])
ylim([-0.02 0.15])
ylabel('dF/F')
xlabel('Off Interval (s)')
title(['Cycle ' num2str(icyc+1)])
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt_quant.pdf']),'-dpdf', '-bestfit')

%find cells responsive to target or at least one target direction
data_dfof_targ = nan(40,nCells,nTrials);
tFramesOff_targ = nan(1,nTrials);
for itrial = 1:nTrials
if ~isnan(cTarget(itrial))
    icyc = nCyc(itrial);
    data_dfof_targ(:,:,itrial) = data_dfof(:,:,icyc,itrial);
    tFramesOff_targ(1,itrial) = tFramesOff(itrial,icyc-1);
end
end
[x_targ y_targ] = ttest(squeeze(nanmean(data_dfof_targ(base_win,:,:),1))',squeeze(nanmean(data_dfof_targ(resp_win,:,:),1))','tail','left','alpha',0.05);

if ndir > 1
h_targ = zeros(ndir,nCells);
p_targ = zeros(ndir,nCells);
for idir = 1:ndir
    ind1 = intersect(find(baseDir == dirs(idir)),find(targetDelta== deltas(1)));
    shift = idir-2;
    if shift<1
        shift = shift+input.baseGratingDirectionStepN;
    end
    ind2 = intersect(find(baseDir == dirs(shift)),find(targetDelta== deltas(2)));
    ind = [ind1 ind2];
    for iCell = 1:nCells
        [h_targ(idir,iCell), p_targ(idir,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof_targ(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
    end
end
else
h_targ = zeros(nDelta,nCells);
p_targ = zeros(nDelta,nCells);
for itarg = 1:nDelta
    ind = find(targetDelta== deltas(itarg));
    for iCell = 1:nCells
        [h_targ(itarg,iCell), p_targ(itarg,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
    end
end
end
good_targ_ind = intersect(good_ind, unique([find(sum(h_targ,1)>0) find(x_targ)]));
good_resp_ind = unique([find(sum(h_targ,1)>0) find(x_targ) good_ind]);


save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']), 'data_dfof', 'max_dir','good_ind', 'good_targ_ind', 'good_resp_ind', 'base_win','resp_win','tt')

figure;
[n n2] = subplotn(length(good_resp_ind));
for iCell = 1:length(good_resp_ind)
    iC = good_resp_ind(iCell);
    subplot(n,n2,iCell)
    plot(tt,data_dfof_base(:,iC,1))
    hold on
    for itarg = 1:nDelta
        ind = find(targetDelta == deltas(itarg));
        plot(tt,nanmean(data_dfof_targ(:,iC,ind),3)-nanmean(nanmean(data_dfof_targ(base_win,iC,ind),1),3));
        hold on;
    end
end
suptitle([mouse ' ' date]) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byCell.pdf']),'-dpdf','-bestfit')

targ_only_ind = unique([find(sum(h_targ,1)>0) find(x_targ)]);
indn = zeros(nDelta, noff);
figure;
for itarg = 1:nDelta
    for ioff = 1:noff
        ind = intersect(find(targetDelta == deltas(itarg)), find(tFramesOff_targ== offs(ioff)));
        indn(itarg,ioff) = length(ind);
        subplot(2,nDelta,itarg)
        plot(tt, nanmean(nanmean(bsxfun(@minus,data_dfof_targ(:,targ_only_ind,ind), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3),2))
        hold on
        subplot(2,nDelta,itarg+nDelta)
        temp = nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,targ_only_ind,ind),1), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3);
        errorbar(offs(ioff)*frameRateHz, nanmean(temp,2), nanstd(temp,[],2)./sqrt(length(targ_only_ind)),'o')
        hold on
    end
    ylim([0 0.4])
    ylabel('dF/F')
    xlabel('Off interval (ms)')
    subplot(2,nDelta,itarg)
    ylim([-0.05 0.4])
    ylabel('dF/F')
    xlabel('Time (ms)')
    title([num2str(deltas(itarg)) ' deg change- ' num2str(indn(itarg,:))])
end
suptitle([mouse ' ' date '- ' num2str(length(targ_only_ind)) ' cells']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byInt_targInd.pdf']),'-dpdf','-bestfit')

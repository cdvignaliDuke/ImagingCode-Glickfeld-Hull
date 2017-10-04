%% Find responsive cells
if length(unique(nCyc)) == 1
    targCyc = unique(nCyc);
else
    targCyc = double(nCyc);
    targCyc(find(FIx)) = NaN;
end
sz = size(data_dfof);
[x1,y1] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
if length(targCyc)==1
    data_dfof_targ = squeeze(data_dfof(:,:,targCyc,:));
else
    data_dfof_targ = nan(sz(1),sz(2),sz(4));
    targInt = nan(1,nTrials);
    for itrial = 1:nTrials
        if ~isnan(targCyc(itrial))
            data_dfof_targ(:,:,itrial) = data_dfof(:,:,targCyc(itrial),itrial);
            targInt(1,itrial) = tFramesOff(itrial,targCyc(itrial)-1);
        end
    end
end
[x2,y2] = ttest(squeeze(mean(data_dfof_targ(base_win,:,:),1))',squeeze(mean(data_dfof_targ(resp_win,:,:),1))','tail','left','alpha',0.05);

        
h1 = zeros(ndir,nCells);
p1 = zeros(ndir,nCells);
h2 = zeros(nDelta,nCells);
p2 = zeros(nDelta,nCells);

data_dfof_dir = nan(sz(1), nCells,ndir);
data_dfof_delta = nan(sz(1), nCells,nDelta);
if ndir>1
    for idir = 1:ndir
        ind = setdiff(find(baseDir == dirs(idir)),ind_con);
        data_dfof_dir(:,:,idir) = nanmean(data_dfof(:,:,1,ind),4);
        for iCell = 1:nCells
            [h1(idir,iCell), p1(idir,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,1,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,1,ind),1)),'tail','left','alpha',0.05./(ndir-1));
        end
    end
else
    data_dfof_dir(:,:,1) = nanmean(data_dfof(:,:,1,:),4);
end
if nDelta>1
    for idelta = 1:nDelta
        ind_delta = find(targetDelta == deltas(idelta));
        ind = setdiff(ind_delta, find(FIx));
        data_dfof_delta(:,:,idelta) = nanmean(data_dfof_targ(:,:,ind),3);
        for iCell = 1:nCells
            [h2(idelta,iCell), p2(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof_targ(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
        end
    end
else
    ind = setdiff(1:nTrials, find(FIx));
    data_dfof_delta(:,:,1) = nanmean(data_dfof_targ(:,:,ind),3);
end

good_ind_temp = unique([find(x1) find(sum(h1,1))]);
good_ind_targ = unique([find(x2) find(sum(h2,1))]);

%find late responding cells and remove
if length(FIx)<1
    tc_all = squeeze(nanmean(mean(bsxfun(@minus,data_dfof(:,:,1,:),mean(data_dfof(base_win,:,1,:),1)),3),4));
    resp_diff = diff(tc_all);
    [max_val max_time] = max(resp_diff(20:31,:),[],1);

    ind1 = find(max_time<3);
    ind2 = find(max_time>9);
    ind3 = setdiff(1:nCells, [ind1 ind2]);
    good_ind = intersect(good_ind_temp,ind3);

    figure;
    subplot(2,2,1)
    plot(tc_all(:,ind1))
    subplot(2,2,2)
    plot(tc_all(:,ind2))
    subplot(2,2,3)
    plot(tc_all(:,ind3))
    subplot(2,2,4)
    plot(tc_all(:,good_ind))
else
    good_ind = good_ind_temp;
end

%% baseline response

tt = (1-20:1+99)*(1000/frameRateHz);

data_dfof_off = cell(noff,noff+1);
data_dfof_off_nminus1 = cell(noff,noff+1);
for ioff = 1:noff
    for icyc = 2:maxCyc-1
        ind = find(tFramesOff(:,icyc) == offs(ioff));
        data_dfof_off{ioff,1} = cat(3, data_dfof_off{ioff,1}, squeeze(data_dfof(:,:,icyc+1,ind)));
        data_dfof_off_nminus1{ioff,1} = cat(3, data_dfof_off{ioff,1}, squeeze(data_dfof(:,:,icyc,ind)));
        for io = 1:noff
            ind_sub = intersect(ind, find(tFramesOff(:,icyc-1) == offs(io)));
            data_dfof_off{ioff,1+io} = cat(3, data_dfof_off{ioff,1+io}, squeeze(data_dfof(:,:,icyc+1,ind_sub)));
            data_dfof_off_nminus1{ioff,1+io} = cat(3, data_dfof_off{ioff,1+io}, squeeze(data_dfof(:,:,icyc,ind_sub)));
        end
    end
end


data_dfof_off_avg = zeros(nCells,noff,noff+1);
data_dfof_off_auroc = zeros(nCells,noff,noff+1);
for ioff  = 1:noff
        temp = data_dfof_off{ioff,1};
        temp_nminus1 = data_dfof_off_nminus1{ioff,1};
%         subplot(2,2,1)
%         plot(tt,nanmean(bsxfun(@minus,temp(:,iC,:),nanmean(temp(base_win,iC,:),1)),3));
%         ylim([-0.1 0.3])
%         hold on
%         subplot(2,2,2)
        temp_avg = bsxfun(@minus,mean(temp(resp_win,iC,:),1),mean(temp(base_win,iC,:),1));
        temp_nminus1_avg = bsxfun(@minus,mean(temp_nminus1(resp_win,iC,:),1),mean(temp_nminus1(base_win,iC,:),1));
%         errorbar(offs(ioff)*(1000/frameRateHz), nanmean(temp_avg,3), nanstd(temp_avg,[],3)./sqrt(size(temp_avg,3)),'o');
%         ylim([-0.05 0.3])
%         hold on
        data_dfof_off_avg(iC,ioff,1) = nanmean(temp_avg,3);
        data_dfof_off_auroc(iC,ioff,1) = roc_gh(temp_nminus1_avg, temp_avg);
%         subplot(2,2,3)
%         plot(offs(ioff)*(1000/frameRateHz),data_dfof_off_auroc(iC,ioff), 'o')
%         hold on
%         title('auROC')
        for io = 1:noff
            temp = data_dfof_off{ioff,1+io};
            temp_nminus1 = data_dfof_off_nminus1{ioff,1+io};
            temp_avg = bsxfun(@minus,mean(temp(resp_win,iC,:),1),mean(temp(base_win,iC,:),1));
            temp_nminus1_avg = bsxfun(@minus,mean(temp_nminus1(resp_win,iC,:),1),mean(temp_nminus1(base_win,iC,:),1));
            data_dfof_off_avg(iC,ioff,1+io) = nanmean(temp_avg,3);
            data_dfof_off_auroc(iC,ioff,1+io) = roc_gh(temp_nminus1_avg, temp_avg);
        end
    end
%     suptitle([mouse ' ' date ' Cell #' num2str(iC)]) 
%     print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_Cell' num2str(iC) '.pdf']),'-dpdf', '-bestfit')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_randIntResp.mat']), 'data_dfof_off_avg', 'data_dfof_off_auroc', 'base_resp', 'good_ind')

%% baseline response figures
tt = (1-20:1+99)*(1000/frameRateHz);
base_dfof = nan(size(data_dfof));
for itrial = 1:nTrials
    n = tCyc(itrial);
    if FIx(itrial)
        n = n-1;
    end
    if n>0
        base_dfof(:,:,1:n,itrial) =data_dfof(:,:,1:n,itrial);
    end
end
base_resp = squeeze(bsxfun(@minus,base_dfof,mean(base_dfof(base_win+3,:,:,:),1)));

figure; 
longHoldTimes = find(celleqel2mat_padded(input.holdTimesMs)>3300);
[n n2] = subplotn(length(good_ind));
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    subplot(n,n2,iCell)
    plot(tt,nanmean(mean(base_resp(:,iC,1,longHoldTimes),2),4));
    hline(0)
    vline(0)
    xlim([tt(1) tt(end)])
    title(['Cell ' num2str(iC)]);
end
suptitle([mouse ' ' date ' - Average anticipation response for anticipation responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgAnticipationResp.pdf']),'-dpdf','-bestfit')

figure;
[n n2] = subplotn(length(good_ind)+1);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    subplot(n,n2,iCell)
    for icyc = 1:maxCyc
        if sum(~isnan(tFramesOff(:,icyc)),1)>20
            plot(tt(1:40),nanmean(mean(base_resp(1:40,iC,icyc,:),2),4));
            hold on
        end
    end
    title(['Cell ' num2str(iC)]);
end
subplot(n,n2,iCell+1)
for icyc = 1:maxCyc
    if sum(~isnan(tFramesOff(:,icyc)),1)>20
        plot(tt(1:40),nanmean(mean(base_resp(1:40,good_ind,icyc,:),2),4));
        hold on
    end
end
title(['All cells']);
suptitle([mouse ' ' date ' - Average response by cycle for anticipation responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgCycleResp.pdf']),'-dpdf','-bestfit')

figure;
[n, n2] = subplotn(length(good_ind)+1);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    subplot(n,n2,iCell)
    for ioff = 1:noff
        base_resp_off = [];
        for icyc = 3:maxCyc
            ind = find(tFramesOff(:,icyc-1) == offs(ioff));
            base_resp_off = cat(4,base_resp_off,base_resp(:,:,icyc,ind));
        end
        plot(tt(1:40),nanmean(mean(base_resp_off(1:40,iC,:,:),2),4))
        hold on
    end
    title(['Cell ' num2str(iC)]);
end
subplot(n,n2,iCell+1)
for ioff = 1:noff
    base_resp_off = [];
    for icyc = 3:maxCyc
        ind = find(tFramesOff(:,icyc-1) == offs(ioff));
        base_resp_off = cat(4,base_resp_off,base_resp(:,:,icyc,ind));
    end
    plot(tt(1:40),nanmean(mean(base_resp_off(1:40,good_ind,:,:),2),4))
    hold on
end
title('All responsive neurons')
suptitle([mouse ' ' date ' - Average response by ISI for cycles 3-' num2str(maxCyc) 'for anticipation responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgISIResp.pdf']),'-dpdf','-bestfit')

%% target response figures
targ_resp = squeeze(bsxfun(@minus,data_dfof_targ,mean(data_dfof_targ(base_win,:,:),1)));

figure; 
[n n2] = subplotn(length(good_ind_targ)+1);
for iCell = 1:length(good_ind_targ)
    iC = good_ind_targ(iCell);
    subplot(n,n2,iCell)
    for idelta = 1:nDelta
        ind = find(tGratingDir == deltas(idelta));
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,iC,ind),2),3));
        hold on
    end
    hline(0)
    vline(0)
    xlim([tt(1) tt(40)])
    title(['Cell ' num2str(iC)]);
    if iCell == 1
        legend(num2str(deltas'))
    end
end
subplot(n,n2,iCell+1)
for idelta = 1:nDelta
    ind = find(tGratingDir == deltas(idelta));
    plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
    hold on
end
hline(0)
vline(0)
xlim([tt(1) tt(40)])
title('All Resp Cells');
suptitle([mouse ' ' date ' - Average response to each target for target responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp.pdf']),'-dpdf','-bestfit')

figure;
[n, n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_delta = find(tGratingDir == deltas(idelta));
    subplot(n,n2,idelta)
    ind_h = intersect(find(SIx),ind_delta);
    plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind_h),2),3));
    hold on
    ind_m = intersect(find(MIx),ind_delta);
    plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind_m),2),3));
    title([num2str(deltas(idelta)) ' deg; H:' num2str(length(ind_h)) ' & M:' num2str(length(ind_m))])
end
suptitle([mouse ' ' date ' - Average response to each target for hits (blue) and misses (red)'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byDelta_byOutcome.pdf']),'-dpdf','-bestfit')

figure;
off_n = zeros(nDelta,noff);
[n, n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_delta = find(tGratingDir == deltas(idelta));
    subplot(n,n2,idelta)
    for ioff = 1:noff
        ind = intersect(ind_delta, find(targInt == offs(ioff)));
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
        hold on
        off_n(idelta,ioff) = length(ind);
    end
    title([num2str(deltas(idelta)) ' deg; n = ' num2str(off_n(idelta,:))])
end
suptitle([mouse ' ' date ' - Average response to each target for all target responsive neurons by ISI'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byDelta_byISI.pdf']),'-dpdf','-bestfit')

figure;
off_n_hit = zeros(nDelta,noff);
[n, n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_delta = find(tGratingDir == deltas(idelta));
    subplot(n,n2,idelta)
    for ioff = 1:noff
        ind = intersect(find(SIx),intersect(ind_delta, find(targInt == offs(ioff))));
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
        hold on
        off_n_hit(idelta,ioff) = length(ind);
    end
    title([num2str(deltas(idelta)) ' deg; n = ' num2str(off_n_hit(idelta,:))])
end
suptitle([mouse ' ' date ' - Average response to Hits for each target for all target responsive neurons by ISI'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byDelta_byISI_hits.pdf']),'-dpdf','-bestfit')

figure;
off_n_miss = zeros(nDelta,noff);
[n, n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_delta = find(tGratingDir == deltas(idelta));
    subplot(n,n2,idelta)
    for ioff = 1:noff
        ind = intersect(find(MIx),intersect(ind_delta, find(targInt == offs(ioff))));
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
        hold on
        off_n_miss(idelta,ioff) = length(ind);
    end
    title([num2str(deltas(idelta)) ' deg; n = ' num2str(off_n_miss(idelta,:))])
end
suptitle([mouse ' ' date ' - Average response to Misses for each target for all target responsive neurons by ISI'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byDelta_byISI_misses.pdf']),'-dpdf','-bestfit')

evenRespMat = cell(1,noff);
for idelta = 1:nDelta
    min_val = min(off_n(idelta,:),[],2);
    ind_delta = find(tGratingDir == deltas(idelta));
    if min_val>0
        for ioff = 1:noff
            ind = intersect(ind_delta, find(targInt == offs(ioff)));
            ind_choose = ind(randsample(length(ind),min_val));
            f(idelta,ioff) = sum(FIx(ind_choose));
            evenRespMat{ioff} = cat(3, evenRespMat{ioff}, targ_resp(:,:,ind_choose));
        end
    end
end

figure; 
[n, n2] = subplotn(length(good_ind_targ)+1); 
for iCell = 1:length(good_ind_targ)
    subplot(n,n2,iCell)
    iC = good_ind_targ(iCell);
    for ioff = 1:noff
        temp = evenRespMat{ioff};
        plot(tt(1:40),nanmean(temp(1:40,iC,:),3))
        hold on
    end
    title(['Cell ' num2str(iC)])
end
subplot(n,n2,iCell+1)
for ioff = 1:noff
    temp = evenRespMat{ioff};
    plot(tt(1:40),nanmean(mean(temp(1:40,good_ind_targ,:),2),3))
    hold on
end
title(['All cells'])
suptitle([mouse ' ' date ' - Average response to each ISI for target responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byISI.pdf']),'-dpdf','-bestfit')

figure;
for ioff = 1:noff
    subplot(1,3,ioff)
    ind_off = find(targInt == offs(ioff));
    for idelta = 1:nDelta
        ind = intersect(ind_off,find(tGratingDir == deltas(idelta)));
        plot(tt(1:40),nanmean(mean(targ_resp(1:40,good_ind_targ,ind),2),3));
        hold on
    end
    xlim([tt(1) tt(40)])
    ylim([-0.02 0.1])
    title([num2str(chop(offs(ioff)*(1000./frameRateHz),3)) ' ms'])
    axis square
end
suptitle([mouse ' ' date ' - Average response to target for each ISI for all target responsive neurons'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTargResp_byISI_byDelta.pdf']),'-dpdf','-bestfit')
            

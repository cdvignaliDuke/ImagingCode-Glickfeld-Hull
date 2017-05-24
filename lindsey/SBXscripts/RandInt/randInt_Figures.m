%% Find responsive cells
if length(unique(nCyc)) == 1
    targCyc = unique(nCyc);
end
[x1,y1] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
[x2,y2] = ttest(squeeze(mean(data_dfof(base_win,:,targCyc,:),1))',squeeze(mean(data_dfof(resp_win,:,targCyc,:),1))','tail','left','alpha',0.05);

h1 = zeros(ndir,nCells);
p1 = zeros(ndir,nCells);
h2 = zeros(nDelta,nCells);
p2 = zeros(nDelta,nCells);

data_dfof_dir = nan(40, nCells,ndir);
data_dfof_delta = nan(40, nCells,nDelta);
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
        ind = find(targetDelta == deltas(idelta));
        data_dfof_delta(:,:,idelta) = nanmean(data_dfof(:,:,targCyc,ind),4);
        for iCell = 1:nCells
            [h2(idelta,iCell), p2(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,targCyc,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,targCyc,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
        end
    end
else
    data_dfof_delta(:,:,1) = nanmean(data_dfof(:,:,targCyc,:),4);
end

good_ind_temp = find(x1);

%find late responding cells and remove
tc_all = squeeze(nanmean(mean(bsxfun(@minus,data_dfof(:,:,1,:),mean(data_dfof(base_win,:,1,:),1)),3),4));
resp_diff = diff(tc_all);
[max_val max_time] = max(resp_diff(20:end,:),[],1);

ind1 = find(max_time<3);
ind2 = find(max_time>9);
ind3 = setdiff(1:nCells, [ind1 ind2]);
good_ind = intersect(good_ind_temp,ind3);

% figure;
% subplot(2,2,1)
% plot(tc_all(:,ind1))
% subplot(2,2,2)
% plot(tc_all(:,ind2))
% subplot(2,2,3)
% plot(tc_all(:,ind3))
% subplot(2,2,4)
% plot(tc_all(:,good_ind))

%%

tt = (1-19:1+20)*(1000/frameRateHz);

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

base_resp = nanmean(bsxfun(@minus,mean(data_dfof(resp_win,:,1,:),1),mean(data_dfof(base_win,:,1,:),1)),4);
data_dfof_off_avg = zeros(nCells,noff,noff+1);
data_dfof_off_auroc = zeros(nCells,noff,noff+1);
for iCell = 1:length(good_ind)
    %figure;
    iC = good_ind(iCell);
%     subplot(2,2,1)
%     plot(tt,nanmean(bsxfun(@minus,data_dfof(:,iC,1,:),nanmean(data_dfof(base_win,iC,1,:),1)),4));
%     hold on
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
end
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_randIntResp.mat']), 'data_dfof_off_avg', 'data_dfof_off_auroc', 'base_resp', 'good_ind')

        
    

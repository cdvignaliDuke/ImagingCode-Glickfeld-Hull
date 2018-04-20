%% Find responsive cells

[x1,y1] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
[x2,y2] = ttest(squeeze(mean(data_dfof(base_win,:,2,:),1))',squeeze(mean(data_dfof(resp_win,:,2,:),1))','tail','left','alpha',0.05);

h1 = zeros(ndir,nCells);
p1 = zeros(ndir,nCells);
h2 = zeros(nDelta,nCells);
p2 = zeros(nDelta,nCells);

data_dfof_dir = nan(size(data_dfof,1), nCells,ndir);
data_dfof_delta = nan(size(data_dfof,1), nCells,nDelta);
if ndir>1
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
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
        data_dfof_delta(:,:,idelta) = nanmean(data_dfof(:,:,2,ind),4);
        for iCell = 1:nCells
            [h2(idelta,iCell), p2(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,2,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,2,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
        end
    end
else
    data_dfof_delta(:,:,1) = nanmean(data_dfof(:,:,2,:),4);
end

good_ind_temp = find(x1+x2+sum(h1,1)+sum(h2,1));

%find max direction for each cell
base_resp = zeros(nCells,ndir);
for iCell = 1:nCells
    for idir = 1:ndir
        ind = find(tGratingDir == dirs(idir));
        base_resp(iCell,idir) = nanmean(mean(data_dfof(resp_win,iCell,1,ind),1)-mean(data_dfof(base_win,iCell,1,ind),1),4);
    end
end

[max_val max_dir] = max(base_resp,[],2);

%find late responding cells and remove
if noff>2
    ind = find(tFramesOff>30);
else
    ind = 1:ntrials;
end
base_tc_all = squeeze(nanmean(bsxfun(@minus,data_dfof(:,:,1,ind),mean(data_dfof(base_win,:,1,ind),1)),4));
base_resp_diff = diff(base_tc_all);
[max_val max_time] = max(base_resp_diff(20:40,:),[],1);

ind1 = find(max_time<base_win(end)-20);
ind2 = find(max_time>resp_win(end)-20);
ind3 = setdiff(1:nCells, [ind1 ind2]);
%ind4 = find(mean(base_tc_all(resp_win,:),1)>0.02);
%good_ind = intersect(intersect(good_ind_temp,ind3),ind4);
good_ind = intersect(good_ind_temp,ind3); %used for V1 PP data 
figure;
subplot(2,2,1)
plot(base_tc_all(:,ind1))
subplot(2,2,2)
plot(base_tc_all(:,ind2))
subplot(2,2,3)
plot(base_tc_all(:,ind3))
subplot(2,2,4)
plot(base_tc_all(:,good_ind))

%% account for blank trials
if length(ind_con)>0
    off_all = [offs; unique(cell2mat(input.tItiWaitFrames))];
    noff_all = length(off_all);
else
    off_all = offs;
    noff_all = noff;
end

%% find fit tau and subtract base response

data_dfof_sub = nan(size(data_dfof,1),nCells,2,nTrials);
ind_long = find(tFramesOff == offs(noff));
data_dfof_sub(:,:,1,:) = data_dfof(:,:,1,:);
data_dfof_sub(:,:,2,ind_long) = data_dfof(:,:,2,ind_long);
R_square = nan(nCells,1);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    fprintf('%d\n',iC)
    tc = nanmean(data_dfof(resp_win(1):end,iC,1,ind_long),4);
    [A, tau, R_square(iC,1)] = singExpDecayFit(tc);
    if isnan(A)
        continue
    end
    fit_resp = A.*exp((1:length(tc))./tau);
    fit_resp_fill = [zeros(1,size(data_dfof,1)-length(fit_resp)) fit_resp];
    for ioff = 1:noff-1
        delay = -offs(ioff)-3;
        fit_shift = circshift(fit_resp_fill',delay);
        ind = find(tFramesOff == off_all(ioff));
        for itrial = 1:length(ind)
            it = ind(itrial);
            base_win_amp = mean(data_dfof(resp_win,iC,1,it),1)-mean(data_dfof(base_win,iC,1,it),1);
            fit_resp_scale = (base_win_amp./fit_resp(1)).*fit_shift;
            data_dfof_sub(:,iC,2,it) = data_dfof(:,iC,2,it)-fit_resp_scale;
        end
    end
end

figure;
[n n2] = subplotn(length(good_ind));
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    subplot(n,n2,iCell)
    ind = setdiff(find(tFramesOff == off_all(1)),ind_long);
    plot(squeeze(nanmean(data_dfof(:,iC,2,ind),4)))
    hold on
    plot(squeeze(nanmean(data_dfof_sub(:,iC,2,ind),4)))
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']), 'data_dfof_sub', 'good_ind', 'off_all', 'noff_all')


%% average responses

temp_resp = squeeze(mean(data_dfof(resp_win,:,:,:),1)-mean(data_dfof(base_win,:,:,:),1));
temp_resp_sub = squeeze(mean(data_dfof_sub(resp_win,:,:,:),1)-mean(data_dfof_sub(base_win,:,:,:),1));

resp_off = zeros(nCells,noff+1);
resp_off_sub = zeros(nCells,noff+1);
resp_off_dir = zeros(nCells,ndir,noff+1);
resp_off_dir_sub = zeros(nCells,ndir,noff+1);
for ioff = 1:noff
    ind = find(tFramesOff == offs(ioff));
    resp_off(:,ioff) = squeeze(nanmean(temp_resp(:,2,ind),3));
    resp_off_sub(:,ioff) = squeeze(nanmean(temp_resp_sub(:,2,ind),3));
    for idir = 1:ndir
        ind_dir = find(baseDir == dirs(idir));
        ind_off = intersect(ind, ind_dir);
        resp_off_dir(:,idir,ioff) = squeeze(nanmean(temp_resp(:,2,ind_off),3));
        resp_off_dir_sub(:,idir,ioff) = squeeze(nanmean(temp_resp_sub(:,2,ind_off),3));
        if ioff == noff
            resp_off_dir(:,idir,noff+1) =  squeeze(nanmean(temp_resp(:,1,ind_dir),3));
            resp_off_dir_sub(:,idir,noff+1) =  squeeze(nanmean(temp_resp_sub(:,1,ind_dir),3));
        end
    end
    if ioff == noff
        resp_off(:,noff+1) =  squeeze(nanmean(temp_resp(:,1,:),3));
        resp_off_sub(:,noff+1) =  squeeze(nanmean(temp_resp_sub(:,1,:),3));
    end
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']), 'resp_off', 'resp_off_sub','resp_off_dir', 'resp_off_dir_sub', 'max_dir', 'good_ind')

%% Figures

tt = (1-20:1+99).*(1000/frameRateHz);
%response by interval- all and max dir- by cell
for iCell = 1:length(good_ind)
    figure;
    iC = good_ind(iCell);
    temp_resp_tc = squeeze(nanmean(data_dfof(resp_win,iC,:,:),1)-nanmean(data_dfof(base_win,iC,:,:),1));
    for idir = 1:ndir
        ind_dir = find(tGratingDir == dirs(idir));
        subplot(3,2,1)
        plot(tt, bsxfun(@minus,nanmean(data_dfof(:,iC,1,ind_dir),4),nanmean(nanmean(data_dfof(base_win,iC,1,ind_dir),4),1)));
        hold on
        subplot(3,2,2)
        errorbar(dirs(idir),nanmean(temp_resp_tc(1,ind_dir),2), nanstd(temp_resp_tc(1,ind_dir),[],2)./sqrt(length(ind_dir)), 'ok');
        hold on
    end
    idir = max_dir(iC);
    subplot(3,2,1)
    title('Base')
    subplot(3,2,2)
    title(['Max dir: ' num2str(dirs(idir))])
    xlabel('Orientation')
    ind_dir = find(tGratingDir == dirs(idir));
    subplot(3,2,3)
    plot(tt, bsxfun(@minus,nanmean(data_dfof(:,iC,1,:),4),nanmean(nanmean(data_dfof(base_win,iC,1,:),4),1)));
    title('All dir')
    hold on
    subplot(3,2,4)
    errorbar(8000, nanmean(temp_resp_tc(1,:),2), nanstd(temp_resp_tc(1,:),[],2)./sqrt(size(data_dfof,4)), 'ob');
    hold on
    xlabel('ISI')
    subplot(3,2,5)
    plot(tt, bsxfun(@minus,nanmean(data_dfof(:,iC,1,ind_dir),4),nanmean(nanmean(data_dfof(base_win,iC,1,ind_dir),4),1)));
    title('Max dir')
    hold on
    subplot(3,2,6)
    errorbar(8000, nanmean(temp_resp_tc(1,ind_dir),2), nanstd(temp_resp_tc(1,ind_dir),[],2)./sqrt(length(ind_dir)), 'ob');
    xlabel('ISI')
    hold on
    for ioff = 1:noff
        ind_off = find(tFramesOff == offs(ioff));
        ind = intersect(ind_dir,ind_off);
        subplot(3,2,3)
        plot(tt, bsxfun(@minus,nanmean(data_dfof(:,iC,2,ind_off),4),nanmean(nanmean(data_dfof(base_win,iC,2,ind_off),4),1)));
        hold on;
        subplot(3,2,4)
        errorbar(offs(ioff)*(1000/frameRateHz), nanmean(temp_resp_tc(2,ind_off),2), nanstd(temp_resp_tc(2,ind_off),[],2)./sqrt(length(ind_off)), 'ob');
        hold on;
        subplot(3,2,5)
        plot(tt, bsxfun(@minus,nanmean(data_dfof(:,iC,2,ind),4),nanmean(nanmean(data_dfof(base_win,iC,2,ind),4),1)));
        hold on;
        subplot(3,2,6)
        errorbar(offs(ioff)*(1000/frameRateHz), nanmean(temp_resp_tc(2,ind),2), nanstd(temp_resp_tc(2,ind),[],2)./sqrt(length(ind)), 'ob');
        hold on;
    end
    suptitle([mouse ' ' date '; Cell #' num2str(iC)])
    %print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_resp_byInt_Cell' num2str(iC) '.pdf']),'-dpdf','-bestfit')
end

%average all cells by interval all directions
tt = (1-20:1+99).*(1000/frameRateHz);
b_tt = (base_win-20).*(1000/frameRateHz);
r_tt = (resp_win-20).*(1000/frameRateHz);
figure;
temp_resp_tc = bsxfun(@minus,data_dfof(:,:,:,:),nanmean(data_dfof(base_win,:,:,:),1));
ind_off = [];
subplot(2,1,1)
plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,1,:),2),4)))
ylabel('dF/F')
hold on
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,:),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,1,:),1),4),[],2)./sqrt(length(good_ind)), 'ob');
ylabel('Average dF/F')
hold on
for ioff = 1:noff
    ind = find(tFramesOff == offs(ioff));
    ind_off = [ind_off length(ind)];
    subplot(2,1,1)
    plot(tt, squeeze(nanmean(nanmean(temp_resp_tc(:,good_ind,2,ind),2),4)))
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind),1),2),4), nanstd(nanmean(nanmean(temp_resp_tc(resp_win,good_ind,2,ind),1),4),[],2)./sqrt(length(good_ind)), 'ob');
    hold on;
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells- ' num2str(ind_off)]) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_allDir.pdf']),'-dpdf','-bestfit')

%average all cells by interval all directions, normalized to base
figure;
norm_resp_tc = bsxfun(@rdivide, temp_resp_tc, nanmean(mean(temp_resp_tc(resp_win,:,1,:),1),4));
ind_off = [];
subplot(2,1,1)
plot(tt, squeeze(nanmean(nanmean(norm_resp_tc(:,good_ind,1,:),2),4)))
ylabel('Normalized dF/F')
hold on
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(nanmean(norm_resp_tc(resp_win,good_ind,1,:),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc(resp_win,good_ind,1,:),1),4),[],2)./sqrt(length(good_ind)), 'ob');
ylabel('Normalized dF/F')
hold on
for ioff = 1:noff
    ind = find(tFramesOff == offs(ioff));
    ind_off = [ind_off length(ind)];
    subplot(2,1,1)
    plot(tt, squeeze(nanmean(nanmean(norm_resp_tc(:,good_ind,2,ind),2),4)))
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(nanmean(norm_resp_tc(resp_win,good_ind,2,ind),1),2),4), nanstd(nanmean(nanmean(norm_resp_tc(resp_win,good_ind,2,ind),1),4),[],2)./sqrt(length(good_ind)), 'ob');
    hold on;
end
ylim([0 1.2])
subplot(2,1,1)
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells- ' num2str(ind_off)]) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_allDir_norm.pdf']),'-dpdf','-bestfit')

%%
%average all cells by interval at maximum direction
data_dfof_maxDir = nan(noff+1,40,nCells);
temp_resp_tc = bsxfun(@minus, data_dfof, nanmean(data_dfof(base_win,:,:,:),1));
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    idir = max_dir(iC);
    ind_dir = find(tGratingDir == dirs(idir));
    data_dfof_maxDir(1,:,iC) =  squeeze(nanmean(temp_resp_tc(:,iC,1,ind_dir),4));
    for ioff = 1:noff
        ind_off = find(tFramesOff == offs(ioff));
        ind = intersect(ind_dir,ind_off);
        data_dfof_maxDir(ioff+1,:,iC) =  squeeze(nanmean(temp_resp_tc(:,iC,2,ind),4));
    end
end

figure;
subplot(2,1,1)
plot(tt, squeeze(nanmean(data_dfof_maxDir(1,:,:),3)))
hold on
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(data_dfof_maxDir(1,resp_win,:),2),3), nanstd(nanmean(data_dfof_maxDir(1,resp_win,:),2),[],3)./sqrt(nCells), 'ob');
hold on
for ioff = 1:noff
    subplot(2,1,1)
    plot(tt, squeeze(nanmean(data_dfof_maxDir(1+ioff,:,:),3)))
    hold on
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(data_dfof_maxDir(1+ioff,resp_win,:),2),3), nanstd(nanmean(data_dfof_maxDir(1+ioff,resp_win,:),2),[],3)./sqrt(nCells), 'ob');
    hold on
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_maxDir.pdf']),'-dpdf','-bestfit')
%%other crap
%average all cells by interval at maximum direction norm to base
figure;
data_dfof_maxDir_norm = bsxfun(@rdivide, data_dfof_maxDir, nanmean(data_dfof_maxDir(1,resp_win,:),2));
subplot(2,1,1)
plot(tt, squeeze(nanmean(data_dfof_maxDir_norm(1,:,:),3)))
hold on
subplot(2,1,2)
errorbar(8000, nanmean(nanmean(data_dfof_maxDir_norm(1,resp_win,:),2),3), nanstd(nanmean(data_dfof_maxDir_norm(1,resp_win,:),2),[],3)./sqrt(nCells), 'ob');
hold on
for ioff = 1:noff
    subplot(2,1,1)
    plot(tt, squeeze(nanmean(data_dfof_maxDir_norm(1+ioff,:,:),3)))
    hold on
    subplot(2,1,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(nanmean(data_dfof_maxDir_norm(1+ioff,resp_win,:),2),3), nanstd(nanmean(data_dfof_maxDir_norm(1+ioff,resp_win,:),2),[],3)./sqrt(nCells), 'ob');
    hold on
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_maxDir_norm.pdf']),'-dpdf','-bestfit')

%scatter of base vs test by interval, average all dirs
[n n2] = subplotn(noff+2);
figure;
temp_resp_tc = bsxfun(@minus, data_dfof, nanmean(data_dfof(base_win,:,:,:),1));
off_resp = zeros(noff+1,nCells);
off_resp_tc = zeros(40,noff+1,nCells);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    off_resp(1,iC) = nanmean(nanmean(temp_resp_tc(resp_win,iC,1,:),1),4);
    off_resp_tc(:,1,iC) = nanmean(temp_resp_tc(:,iC,1,ind_dir),4);
    for ioff = 1:noff
        subplot(n,n2,ioff)
        ind = find(tFramesOff == offs(ioff));
        plot(off_resp(1,iC),nanmean(nanmean(temp_resp_tc(resp_win,iC,2,ind),1),4),'ob')
        hold on;
        off_resp(1+ioff,iC) = nanmean(nanmean(temp_resp_tc(resp_win,iC,2,ind),1),4);
        off_resp_tc(:,1+ioff,iC) = nanmean(temp_resp_tc(:,iC,2,ind),4);
    end
end
for ioff = 1:noff
    subplot(n,n2,ioff)
    ylim([-0.01 0.15])
    xlim([-0.01 0.15])
    xlabel('Base resp dF/F')
    ylabel('Test resp dF/F')
    refline(1,0)
    title(num2str(offs(ioff)*(1000/frameRateHz)))
    errorbarxy(mean(off_resp(1,:),2),mean(off_resp(ioff+1,:),2),std(off_resp(1,:),[],2)./sqrt(length(good_ind)),std(off_resp(ioff+1,:),[],2)./sqrt(length(good_ind)), {'ok', 'k', 'k'})
    hold on
end
off_resp(find(off_resp<0)) = 0;
off_resp_norm = bsxfun(@rdivide,off_resp,off_resp(1,:));
for ioff = 1:noff
    subplot(n,n2,6)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,:),2), nanstd(off_resp_norm(ioff+1,:),[],2)./sqrt(length(good_ind)),'ob')
    hold on
    subplot(n,n2,7)
    plot(tt,nanmean(off_resp_tc(:,ioff+1,:),3))
    hold on
end
subplot(n,n2,6)
ylim([0 1.5])
xlabel('ISI')
subplot(n,n2,7)
xlabel('Time (ms)')
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_scatter.pdf']),'-dpdf','-bestfit')

%heatmap of average all dirs, by interval by cell
figure;
subplot(1,6,1)
imagesc(nanmean(temp_resp_tc(:,:,1,:),4)')
colormap(redblue)
clim([-0.15 0.15])
title('Base')
for ioff = 1:noff
    subplot(1,6,ioff+1)
    ind = find(tFramesOff == offs(ioff));
    imagesc(nanmean(temp_resp_tc(:,:,2,ind),4)')
    colormap(redblue)
    clim([-0.15 0.15])
    title(num2str(offs(ioff)*(1000/frameRateHz)))
end
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_heatmap.pdf']),'-dpdf','-bestfit')

%scatter of base vs test by interval, only for max dir
[n n2] = subplotn(noff+2);
figure;
temp_resp_tc = bsxfun(@minus, data_dfof, nanmean(data_dfof(base_win,:,:,:),1));
off_resp = zeros(noff+1,nCells);
off_resp_tc = zeros(40,noff+1,nCells);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    ind_dir = find(tGratingDir == dirs(max_dir(iC,:)));
    off_resp(1,iC) = nanmean(nanmean(temp_resp_tc(resp_win,iC,1,ind_dir),1),4);
    off_resp_tc(:,1,iC) = nanmean(temp_resp_tc(:,iC,1,ind_dir),4);
    for ioff = 1:noff
        subplot(n,n2,ioff)
        ind = intersect(ind_dir, find(tFramesOff == offs(ioff)));
        plot(off_resp(1,iC),nanmean(nanmean(temp_resp_tc(resp_win,iC,2,ind),1),4),'ob')
        hold on;
        off_resp(1+ioff,iC) = nanmean(nanmean(temp_resp_tc(resp_win,iC,2,ind),1),4);
        off_resp_tc(:,1+ioff,iC) = nanmean(temp_resp_tc(:,iC,2,ind),4);
    end
end
for ioff = 1:noff
    subplot(n,n2,ioff)
    ylim([-0.05 0.25])
    xlim([-0.05 0.25])
    xlabel('Base resp dF/F')
    ylabel('Test resp dF/F')
    refline(1,0)
    title(num2str(offs(ioff)*(1000/frameRateHz)))
    errorbarxy(mean(off_resp(1,:),2),mean(off_resp(ioff+1,:),2),std(off_resp(1,:),[],2)./sqrt(length(good_ind)),std(off_resp(ioff+1,:),[],2)./sqrt(length(good_ind)), {'ok', 'k', 'k'})
    hold on
end
off_resp(find(off_resp<0)) = 0;
off_resp_norm = bsxfun(@rdivide,off_resp,off_resp(1,:));
off_resp_tc_norm = zeros(40,nCells,noff+1);
for ioff = 1:noff+1
    off_resp_tc_norm(:,:,ioff) = bsxfun(@rdivide,squeeze(off_resp_tc(:,ioff,:)),off_resp(1,:));
end
off_resp_tc_norm = permute(off_resp_tc_norm,[1 3 2]);
for ioff = 1:noff
    subplot(n,n2,6)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,:),2), nanstd(off_resp_norm(ioff+1,:),[],2)./sqrt(length(good_ind)),'ob')
    hold on
    subplot(n,n2,7)
    plot(tt,nanmean(off_resp_tc_norm(:,ioff+1,:),3))
    hold on
end
subplot(n,n2,6)
ylim([0 1.5])
xlabel('ISI')
subplot(n,n2,7)
xlabel('Time (ms)')
suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_maxDir_scatter.pdf']),'-dpdf','-bestfit')

%divide cells into facilitating and depressing- using ttest
[h_f p_f] = ttest(mean(temp_resp_tc(resp_win,:,1,:),1), mean(temp_resp_tc(resp_win,:,2,:),1),'tail','left', 'dim',4);
[h_d p_d] = ttest(mean(temp_resp_tc(resp_win,:,1,:),1), mean(temp_resp_tc(resp_win,:,2,:),1),'tail','right','dim',4);

ind_f = intersect(good_ind, find(h_f));
ind_d = intersect(good_ind, find(h_d));

figure;
subplot(2,2,2)
for ioff = 1:noff
    subplot(2,2,1)
    plot(tt,nanmean(off_resp_tc_norm(:,ioff+1,ind_f),3))    
    hold on
    subplot(2,2,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,ind_f),2), nanstd(off_resp_norm(ioff+1,ind_f),[],2)./sqrt(length(ind_f)),'ob')
    hold on
end
subplot(2,2,1)
ylim([-.5 1.5])
title(['Faclitating cells- n = ' num2str(length(ind_f))])
subplot(2,2,2)
ylim([0 1.5])
for ioff = 1:noff
    subplot(2,2,3)
    plot(tt,nanmean(off_resp_tc_norm(:,ioff+1,ind_d),3))    
    hold on
    subplot(2,2,4)
	errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,ind_d),2), nanstd(off_resp_norm(ioff+1,ind_d),[],2)./sqrt(length(ind_d)),'ob')
    hold on
end
subplot(2,2,3)
ylim([-.5 1.5])
title(['Depressing cells- n = ' num2str(length(ind_d))])
subplot(2,2,4)
ylim([0 1.5])
suptitle([mouse ' ' date])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_FDttest.pdf']),'-dpdf','-bestfit')

%divide cells into facilitating
ind_f = intersect(good_ind, find(off_resp_norm(1,:)<off_resp_norm(2,:)));
ind_d = intersect(good_ind, find(off_resp_norm(1,:)>off_resp_norm(2,:)));

figure;
subplot(2,2,2)
for ioff = 1:noff
    subplot(2,2,1)
    plot(tt,nanmean(off_resp_tc_norm(:,ioff+1,ind_f),3))    
    hold on
    subplot(2,2,2)
    errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,ind_f),2), nanstd(off_resp_norm(ioff+1,ind_f),[],2)./sqrt(length(ind_f)),'ob')
    hold on
end
subplot(2,2,1)
ylim([-.5 1.5])
title(['Faclitating cells- n = ' num2str(length(ind_f))])
subplot(2,2,2)
ylim([0 1.5])
for ioff = 1:noff
    subplot(2,2,3)
    plot(tt,nanmean(off_resp_tc_norm(:,ioff+1,ind_d),3))    
    hold on
    subplot(2,2,4)
	errorbar(offs(ioff)*(1000/frameRateHz), nanmean(off_resp_norm(ioff+1,ind_d),2), nanstd(off_resp_norm(ioff+1,ind_d),[],2)./sqrt(length(ind_d)),'ob')
    hold on
end
subplot(2,2,3)
ylim([-.5 1.5])
title(['Depressing cells- n = ' num2str(length(ind_d))])
subplot(2,2,4)
ylim([0 1.5])
suptitle([mouse ' ' date])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_FD.pdf']),'-dpdf','-bestfit')


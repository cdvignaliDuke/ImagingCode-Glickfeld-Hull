%% Find responsive cells

if input.doRandCon
    baseCon = nan(maxCyc,nTrials);
    for itrial = 1:nTrials
        baseCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
    end
end

targCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(targCon);
ncon = length(cons);

ind_base = find(baseCon(1,:)==max(cons));
[h1,p1] = ttest(squeeze(mean(data_dfof(base_win,:,1,ind_base),1))',squeeze(mean(data_dfof(resp_win,:,1,ind_base),1))','tail','left','alpha',0.05);
[h,p] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
good_ind_temp = find(h1+h);

ht1 = zeros(nDelta,nCells);
pt1 = zeros(nDelta,nCells);


targCyc = unique(celleqel2mat_padded(input.tCyclesOn))+1;
if length(targCyc==1)
    if nDelta>1
        for idelta = 1:nDelta
            ind = find(targetDelta == deltas(idelta));
            for iCell = 1:nCells
                [ht1(idelta,iCell), pt1(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,targCyc,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,targCyc,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
            end
        end
    end
    [ht,pt] = ttest(squeeze(mean(data_dfof(base_win,:,targCyc,:),1))',squeeze(mean(data_dfof(resp_win,:,targCyc,:),1))','tail','left','alpha',0.05);
end
good_ind_targ = find(sum([ht;ht1],1));

% ind = find(tFramesOff(:,1)==offs(end));
% tc_all = squeeze(nanmean(bsxfun(@minus,data_dfof(:,:,1,ind),mean(data_dfof(base_win,:,1,ind),1)),4));
% resp_diff = diff(tc_all);
% [max_val max_time] = max(resp_diff(20:50,:),[],1);
% 
% ind1 = find(max_time<base_win(end)-20);
% ind2 = find(max_time>resp_win(end)-20);
% ind3 = setdiff(1:nCells, [ind1 ind2]);
% good_ind_base = intersect(good_ind_temp,ind3);

good_ind_base = good_ind_temp;
good_ind = unique([good_ind_base good_ind_targ]);
good_int = intersect(good_ind_base, good_ind_targ);
good_sub = find(ismember(good_ind_base,good_int));

save(fullfile(LG_base,'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respSig.mat']), 'p','h','p1','h1','pt','ht','pt1','ht1','good_ind','good_ind_base','good_ind_targ')
%% FS rand con
tt= [-21:88].*(1000/frameRateHz);
ttAxisTick = find(ismember(floor(tt),0:500:3000));
ttAxisLabel = floor(tt(ttAxisTick));

[avg out] = sort(nanmean(nanmean(data_notarg_dfof(70:100,good_ind,:),1),3),2,'descend');
data_notarg_dfof_sort = data_notarg_dfof(:,good_ind(out),:);
figure; h = imagesc(nanmean(data_notarg_dfof_sort(1:100,:,:),3)');
figXAxis([],'time (ms)',[],ttAxisTick,ttAxisLabel)
figYAxis([],'cell #',[])
figAxForm([])
colormap('bluered')
colorbar
clim([-0.1 0.1])
title([mouse ' ' date ' ' run_str ': FS rand int; avg all trials'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTCHeatmap.pdf']),'-dpdf','-fillpage')


offs = unique(tFramesOff);
noff = length(offs);

%% sort by n
resp = cell(ncon,noff);
for icon = 1:ncon
    for ioff = 1:noff
        resp{icon,ioff} = [];
        for itrial = 1:nTrials
            ind = intersect(find(baseCon(3:end,itrial)==cons(icon)),find(tFramesOff(itrial,2:end-1)==offs(ioff)));
            if length(ind)>0
                resp{icon,ioff} = cat(3,resp{icon,ioff},data_dfof(:,good_ind_base,ind+2,itrial));
            end
        end
    end
end
resp_off = cell(1,noff);
for ioff = 1:noff
    resp_off{ioff} = [];
    for icon = 1:ncon
        resp_off{ioff} = cat(3,resp_off{ioff},resp{icon,ioff});
    end
end

resp_con = cell(1,ncon);
for icon = 1:icon
    resp_con{icon} = [];
    for ioff = 1:noff
        resp_con{icon} = cat(3,resp_con{icon},resp{icon,ioff});
    end
end

%% amplitude distributions

resp_con_amp = cell(ncon,nDelta+1);
cell_con = zeros(ncon,nDelta+1);
cell_ori = zeros(ncon,nDelta+1);
for icon = 1:ncon
    resp_con_amp{icon,1} = squeeze(mean(resp_con{1,icon}(resp_win,good_sub,:),1)-mean(resp_con{1,icon}(base_win,good_sub,:),1));
    cell_con(icon,1) = cons(icon);
    cell_ori(icon,1) = 0;
    for itarg = 1:nDelta
        ind = intersect(find(tGratingDir == deltas(itarg)), find(targCon == cons(icon)));
        resp_con_amp{icon,itarg+1} = squeeze(mean(data_dfof(resp_win,good_int,end,ind),1)-mean(data_dfof(base_win,good_int,end,ind),1));
        cell_con(icon,itarg+1) = cons(icon);
        cell_ori(icon,itarg+1) = deltas(itarg);
    end
end

save(fullfile(LG_base,'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respDist.mat']),'resp_con_amp','cell_con','cell_ori','offs')

%% cons and offs figures
figure;
subplot(1,2,1)
for icon = 1:ncon
    plot(1:30,nanmean(mean(bsxfun(@minus,resp_con{icon}(11:40,:,:),mean(resp_con{icon}(base_win,:,:),1)),2),3))
    hold on
end
ylim([-0.01 0.04])
title('All Offs')
legend(num2str(cons'))
subplot(1,2,2)
for ioff = 1:noff
    plot(1:30,nanmean(mean(bsxfun(@minus,resp_off{ioff}(11:40,:,:),mean(resp_off{ioff}(base_win,:,:),1)),2),3))
    hold on
end
ylim([-0.01 0.04])
title('All Cons')
legend(num2str(offs.*(1000/frameRateHz)))
suptitle([mouse ' ' date '- Stims 3-' num2str(size(baseCon,1))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgBase_allConsAllOffs.pdf']),'-dpdf','-bestfit')

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))

figure;
resp_con_avg = zeros(ncon,length(good_ind_base));
sz = size(mask_cell);
mask_cell_reshape = reshape(mask_cell,[sz(1)*sz(2) 1]);
mask_resp_con = zeros(sz(1)*sz(2),ncon);
mask_thresh_con = zeros(sz(1)*sz(2),ncon);
[n n2] = subplotn(nCells);
for iCell = 1:length(good_ind_base)
    subplot(n,n2,iCell)
    for icon = 1:ncon
        plot(1:30,nanmean(mean(bsxfun(@minus,resp_con{icon}(11:40,iCell,:),mean(resp_con{icon}(base_win,iCell,:),1)),2),3))
        hold on
        resp_con_avg(icon,iCell) = nanmean(mean(bsxfun(@minus,mean(resp_con{icon}(resp_win,iCell,:)),mean(resp_con{icon}(base_win,iCell,:),1)),2),3);
        
        if resp_con_avg(icon,iCell) > 0
            if ttest(mean(resp_con{icon}(resp_win,iCell,:)),mean(resp_con{icon}(base_win,iCell,:),1),'tail','right')
                mask_resp_con(find(mask_cell_reshape==good_ind_base(iCell)),icon) = resp_con_avg(icon,iCell);
                mask_thresh_con(find(mask_cell_reshape==good_ind_base(iCell)),icon) = 1;
            end
        end
    end
end
mask_resp_con = reshape(mask_resp_con,[sz(1) sz(2) ncon]);
mask_thresh_con = reshape(mask_thresh_con,[sz(1) sz(2) ncon]);

figure;
mask_resp_off = zeros(sz(1)*sz(2),noff);
mask_thresh_off = zeros(sz(1)*sz(2),noff);
resp_off_avg = zeros(noff,length(good_ind_base));
[n n2] = subplotn(nCells);
for iCell = 1:length(good_ind_base)
    subplot(n,n2,iCell)
    for ioff = 1:noff
        plot(1:30,nanmean(mean(bsxfun(@minus,resp_off{ioff}(11:40,iCell,:),mean(resp_off{ioff}(base_win,iCell,:),1)),2),3))
        hold on
        resp_off_avg(ioff,iCell) = nanmean(mean(bsxfun(@minus,mean(resp_off{ioff}(resp_win,iCell,:)),mean(resp_off{ioff}(base_win,iCell,:),1)),2),3);
        if resp_off_avg(ioff,iCell) > 0
            if ttest(mean(resp_off{ioff}(resp_win,iCell,:)),mean(resp_off{ioff}(base_win,iCell,:),1),'tail','right')
                mask_thresh_off(find(mask_cell_reshape==good_ind_base(iCell)),ioff) = 1;
                mask_resp_off(find(mask_cell_reshape==good_ind_base(iCell)),ioff) = resp_off_avg(ioff,iCell);
            end
        end
    end
end
mask_resp_off = reshape(mask_resp_off,[sz(1) sz(2) noff]);
mask_thresh_off = reshape(mask_thresh_off,[sz(1) sz(2) noff]);

figure;
start = 0;
for icon = 1:ncon
    subplot(ncon,2,icon+start)
    imagesc(mask_resp_con(:,:,icon))
    clim([0 0.15])
    title(['Con = ' num2str(cons(icon))])
    start = start+1;
end
start = 1;
for ioff = 1:noff
    subplot(ncon,2,ioff+start)
    imagesc(mask_resp_off(:,:,ioff))
    clim([0 0.15])
    title(['ISI = ' num2str(chop(offs(ioff).*(1000/frameRateHz),3))])
    start = start+1;
end
suptitle([mouse ' ' date '- Response amplitude'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respFOV_dfof.pdf']),'-dpdf','-bestfit')

figure;
start = 0;
for icon = 1:ncon
    subplot(ncon,2,icon+start)
    imagesc(mask_thresh_con(:,:,icon))
    title(['Con = ' num2str(cons(icon)) ' - ' num2str(max(max(bwlabel(mask_thresh_con(:,:,icon)),[],1),[],2)) ' cells'])
    start = start+1;
end
start = 1;
for ioff = 1:noff
    subplot(ncon,2,ioff+start)
    imagesc(mask_thresh_off(:,:,ioff))
    title(['ISI = ' num2str(chop(offs(ioff).*(1000/frameRateHz),3)) ' - ' num2str(max(max(bwlabel(mask_thresh_off(:,:,ioff)),[],1),[],2)) ' cells'])
    start = start+1;
end
suptitle([mouse ' ' date '- Responsive by ttest'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respFOV_thresh.pdf']),'-dpdf','-bestfit')
    
figure;
for ioff = 1:noff
    subplot(1,noff,ioff)
    for icon = 1:ncon
        plot(1:30, nanmean(mean(bsxfun(@minus,resp{icon,ioff}(11:40,:,:),mean(resp{icon,ioff}(base_win,:,:),1)),2),3))
        hold on
    end
    ylim([-0.01 0.04])
    title(['ISI = ' num2str(offs(ioff).*(1000/frameRateHz)) ' ms'])
    if ioff == 1 
        legend(num2str(cons'),'Location','northeast')
    end
end
suptitle([mouse ' ' date '- Stims 3-' num2str(size(baseCon,1))])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgBase_byConsByOffs.pdf']),'-dpdf','-bestfit')

figure;
for icon = 1:ncon
    ind = find(baseCon(1,:) == cons(icon));
    plot(1:30, squeeze(nanmean(mean(data_dfof(11:40,good_ind_base,1,ind),2),4)))
    hold on
end
legend(num2str(cons'))
title([mouse ' ' date '- First stim'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgFirstBase_byCons.pdf']),'-dpdf','-bestfit')

figure;
for ioff = 1:noff
    ind = find(tFramesOff(:,1)==offs(ioff));
    plot(1:40,squeeze(nanmean(mean(data_dfof(1:40,good_ind_base,2,ind)-mean(data_dfof(base_win,good_ind_base,2,ind),1),2),4)));
    hold on
end
legend([num2str(offs.*(1000/frameRateHz))])
title([mouse ' ' date '- First stim'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgSecondBase_byOffs.pdf']),'-dpdf','-bestfit')

%% sort by n-1

resp_n1 = cell(ncon,noff);
for icon = 1:ncon
    for ioff = 1:noff
        resp_n1{icon,ioff} = [];
        for itrial = 1:nTrials
            ind = intersect(find(baseCon(3:end-1,itrial)==cons(icon)),find(tFramesOff(itrial,2:end-2)==offs(ioff)));
            if length(ind)>0
                resp_n1{icon,ioff} = cat(3,resp_n1{icon,ioff},data_dfof(:,good_ind_base,ind+3,itrial));
            end
        end
    end
end
resp_off_n1 = cell(1,noff);
for ioff = 1:noff
    resp_off_n1{ioff} = [];
    for icon = 1:ncon
        resp_off_n1{ioff} = cat(3,resp_off_n1{ioff},resp_n1{icon,ioff});
    end
end

resp_con_n1 = cell(1,ncon);
for icon = 1:icon
    resp_con_n1{icon} = [];
    for ioff = 1:noff
        resp_con_n1{icon} = cat(3,resp_con_n1{icon},resp_n1{icon,ioff});
    end
end

figure;
subplot(1,2,1)
for icon = 1:ncon
    plot(1:30,nanmean(mean(bsxfun(@minus,resp_con_n1{icon}(11:40,:,:),mean(resp_con_n1{icon}(base_win,:,:),1)),2),3))
    hold on
end
ylim([-0.01 0.04])
title('All Offs')
legend(num2str(cons'))
subplot(1,2,2)
for ioff = 1:noff
    plot(1:30,nanmean(mean(bsxfun(@minus,resp_off_n1{ioff}(11:40,:,:),mean(resp_off_n1{ioff}(base_win,:,:),1)),2),3))
    hold on
end
ylim([-0.01 0.04])
title('All Cons')
legend(num2str(offs.*(1000/frameRateHz)))
suptitle([mouse ' ' date '- Stims 3-' num2str(size(baseCon,1)) '; N-1'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgBase_allConsAllOffs_n1.pdf']),'-dpdf','-bestfit')

figure;
for ioff = 1:noff
    subplot(1,noff,ioff)
    for icon = 1:ncon
        plot(1:30, nanmean(mean(bsxfun(@minus,resp_n1{icon,ioff}(11:40,:,:),mean(resp_n1{icon,ioff}(base_win,:,:),1)),2),3))
        hold on
    end
    ylim([-0.01 0.04])
    title(['ISI = ' num2str(offs(ioff).*(1000/frameRateHz)) ' ms'])
    if ioff == 1 
        legend(num2str(cons'),'Location','northeast')
    end
end
suptitle([mouse ' ' date '- Stims 3-' num2str(size(baseCon,1)) '; N-1'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgBase_byConsByOffs_n1.pdf']),'-dpdf','-bestfit')

%for each cell
figure;
resp_con_avg_n1 = zeros(ncon,length(good_ind_base));
sz = size(mask_cell);
mask_cell_reshape = reshape(mask_cell,[sz(1)*sz(2) 1]);
mask_resp_con_n1 = zeros(sz(1)*sz(2),ncon);
[n n2] = subplotn(length(good_ind_base));
for iCell = 1:length(good_ind_base)
    subplot(n,n2,iCell)
    for icon = 1:ncon
        plot(1:30,nanmean(mean(bsxfun(@minus,resp_con_n1{icon}(11:40,iCell,:),mean(resp_con_n1{icon}(base_win,iCell,:),1)),2),3))
        hold on
        resp_con_avg_n1(icon,iCell) = nanmean(mean(bsxfun(@minus,mean(resp_con_n1{icon}(resp_win,iCell,:)),mean(resp_con_n1{icon}(base_win,iCell,:),1)),2),3);
        if resp_con_avg_n1(icon,iCell) > 0
            if ttest(mean(resp_con_n1{icon}(resp_win,iCell,:)),mean(resp_con_n1{icon}(base_win,iCell,:),1),'tail','right')
                mask_resp_con_n1(find(mask_cell_reshape==good_ind_base(iCell)),icon) = resp_con_avg_n1(icon,iCell);
            end
        end
    end
end
mask_resp_con_n1 = reshape(mask_resp_con_n1,[sz(1) sz(2) ncon]);

figure;
mask_resp_off_n1 = zeros(sz(1)*sz(2),noff);
resp_off_avg_n1 = zeros(noff,length(good_ind_base));
[n n2] = subplotn(length(good_ind_base));
for iCell = 1:length(good_ind_base)
    subplot(n,n2,iCell)
    for ioff = 1:noff
        plot(1:30,nanmean(mean(bsxfun(@minus,resp_off_n1{ioff}(11:40,iCell,:),mean(resp_off_n1{ioff}(base_win,iCell,:),1)),2),3))
        hold on
        resp_off_avg_n1(ioff,iCell) = nanmean(mean(bsxfun(@minus,mean(resp_off_n1{ioff}(resp_win,iCell,:)),mean(resp_off_n1{ioff}(base_win,iCell,:),1)),2),3);
        if resp_off_avg_n1(ioff,iCell) > 0
            if ttest(mean(resp_off_n1{ioff}(resp_win,iCell,:)),mean(resp_off_n1{ioff}(base_win,iCell,:),1),'tail','right')
                mask_resp_off_n1(find(mask_cell_reshape==good_ind_base(iCell)),ioff) = resp_off_avg_n1(ioff,iCell);
            end
        end
    end
end
mask_resp_off_n1 = reshape(mask_resp_off_n1,[sz(1) sz(2) noff]);

figure;
start = 0;
for icon = 1:ncon
    subplot(ncon,2,icon+start)
    imagesc(mask_resp_con_n1(:,:,icon))
    clim([0 0.15])
    title(['N-1 Con = ' num2str(cons(icon))])
    start = start+1;
end
start = 1;
for ioff = 1:noff
    subplot(ncon,2,ioff+start)
    imagesc(mask_resp_off_n1(:,:,ioff))
    clim([0 0.15])
    title(['N-1 ISI = ' num2str(chop(offs(ioff).*(1000/frameRateHz),3))])
    start = start+1;
end
suptitle([mouse ' ' date '- Response amplitude'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respFOV_dfof_n1.pdf']),'-dpdf','-bestfit')


%% targets
figure
for itarg = 1:nDelta
    ind = find(tGratingDir == deltas(itarg));
    subplot(2,nDelta,itarg)
    for icon = 1:ncon
        ind2 = intersect(ind, find(targCon == cons(icon)));
        indc(itarg,icon) = length(ind2);
        
        plot(1:30, nanmean(mean(bsxfun(@minus,data_dfof(11:40,good_ind_targ,end,ind2),mean(data_dfof(base_win,good_ind_targ,end,ind2),1)),2),4))
        hold on
        ylim([-0.01 0.06])
    end
    if itarg == 1
        legend(num2str(cons'),'Location','northeast')
    end
    title([num2str(deltas(itarg)) ' deg'])
    subplot(2,nDelta,itarg+nDelta)
    for ioff = 1:noff
        ind2 = intersect(ind, find(tFramesOff(:,end)==offs(ioff)));
        indo(itarg,ioff) = length(ind2);
        plot(1:30, nanmean(mean(bsxfun(@minus,data_dfof(11:40,good_ind_targ,end,ind2),mean(data_dfof(base_win,good_ind_targ,end,ind2),1)),2),4))
        hold on
        ylim([-0.01 0.06])
    end
    if itarg == 1
        legend(num2str(offs.*(1000/frameRateHz)),'Location','northeast')
    end
end
suptitle([mouse ' ' date '- Target by ISI or Contrast'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_targByConOrISI.pdf']),'-dpdf','-bestfit')


figure;
for itarg = 1:nDelta
    for ioff = 1:noff
        ind = intersect(find(tGratingDir == deltas(itarg)), find(tFramesOff(:,end)==offs(ioff)));
        subplot(nDelta,noff,itarg+(nDelta*(ioff-1)))
        for icon = 1:ncon
            ind2 = intersect(ind, find(targCon == cons(icon)));
            n(ioff,icon) = length(ind2);
            plot(1:30, nanmean(mean(bsxfun(@minus,data_dfof(11:40,good_ind_targ,end,ind2),mean(data_dfof(base_win,good_ind_targ,end,ind2),1)),2),4))
            hold on
            ylim([-0.01 0.06])
        end
        legend([num2str(cons') repmat('; n= ',[3,1]) num2str(n(ioff,:)')],'Location','northwest')
        title(['Target = ' num2str(deltas(itarg)) '; ISI = ' num2str(chop(offs(ioff).*(1000/frameRateHz),2)) ' ms'])
    end
end
suptitle([mouse ' ' date '- Target by ISI and Contrast'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_targByConAndISI.pdf']),'-dpdf','-bestfit')





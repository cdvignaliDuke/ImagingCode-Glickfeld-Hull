%% current datasets

mouse_mat = strvcat('i674', 'i689', 'i696','i684','i711','i712','i574');
date_mat = strvcat('170324', '170323', '170323','170327','170503','170503','170510');
ImgFolder = strvcat('002', '003');
nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

%%
nexp = size(mouse_mat,1);
ind_min = cell(1,nexp);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    
    ind{iexp} = zeros(nDelta,noff);
    for idelta = 1:nDelta
        for ioff = 1:noff
            ind{iexp}(idelta,ioff) = length(intersect(find(tFramesOff(:,end) == offs(ioff)), find(targetDelta == deltas(idelta))));
        end
    end
    [ind_min{iexp} i] = min(min(ind{iexp},[],1),[],2);
end

%% collect datasets
nexp = size(mouse_mat,1);
max_dir_all = [];
pref_ori_all = [];
theta_90_all = [];
max_dir_n_all = [];
fit_all = [];
OSI_all = [];
OSI_k_all = [];
HWHM_all = [];
delta_resp_all = [];
roc_resp_all = [];
pref_ori_sub_all = [];
fit_sub_all = [];
OSI_sub_all = [];
OSI_k_sub_all = [];
HWHM_sub_all = [];
delta_resp_sub_all = [];
roc_resp_sub_all = [];
nCells = [];
ppResp_all = cell(1,nexp);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_deltaResp.mat']))
    max_dir_all = [max_dir_all; max_dir];
    pref_ori_all = [pref_ori_all; pref_ori];
    theta_90_all = [theta_90_all; theta_90'];
    max_dir_n_all = [max_dir_n_all; max_dir_n'];
    OSI_all = [OSI_all; OSI];
    OSI_k_all = [OSI_k_all; OSI_k];
    HWHM_all = [HWHM_all; HWHM];
    delta_resp_all = cat(1,delta_resp_all, delta_resp(:,:,:,1));
    fit_all = cat(2, fit_all, y_fit(:,:,:,1));
    pref_ori_sub_all = [pref_ori_sub_all; pref_ori_sub];
    OSI_sub_all = [OSI_sub_all; OSI_sub];
    OSI_k_sub_all = [OSI_k_sub_all; OSI_k_sub];
    HWHM_sub_all = [HWHM_sub_all; HWHM_sub];
    delta_resp_sub_all = cat(1,delta_resp_sub_all, delta_resp_sub(:,:,:,1));
    fit_sub_all = cat(2, fit_sub_all, y_fit_sub(:,:,:,1));
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_roc180v23.mat']))
    roc_resp_all = cat(1,roc_resp_all, roc_resp);
    roc_resp_sub_all = cat(1,roc_resp_sub_all, roc_resp_sub);
    ppResp_all{iexp} = ppResp;
end
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
if length(ind_con)>0
    off_all = [offs; unique(cell2mat(input.tItiWaitFrames))];
    noff_all = length(off_all);
else
    off_all = offs;
    noff_all = noff;
end
good_ind_theta = find(theta_90_all<= 22.5);
good_ind_dir = find(max_dir_n_all>= 500);

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'ppDiffOri_summary.mat'),'mouse_mat', 'date_mat', 'max_dir_all', 'pref_ori_all', 'theta_90_all', 'OSI_all', 'delta_resp_all', 'fit_all', 'good_ind_theta');

%% summary of changes in tuning
figure;
col_mat = strvcat('b', 'r', 'y');
start = 1;
rep = 1;
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    if start>16
        suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Cells ' num2str(iCell-16) ' to ' num2str(iCell-1) '- all fits'])
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['goodFits_summary' num2str(rep) '.pdf']),'-dpdf','-fillpage')
        rep = 1+rep;
        figure;
        start = 1;
    end
    subplot(4, 4, start)
    for ioff = 1:noff_all
        plot(squeeze(fit_all(:,iC,ioff)),col_mat(ioff))
        hold on
        scatter(deltas, delta_resp_all(iC,:,ioff),['o' col_mat(ioff)])
    end
    start = start+1;
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Cells ' num2str(iCell-start+2) ' to ' num2str(iCell) '- all fits'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['goodFits_summary' num2str(rep) '.pdf']),'-dpdf','-fillpage')   

%% tuning figures
delta_resp_norm_all = bsxfun(@rdivide, delta_resp_all, max(delta_resp_all(:,:,3),[],2));
figure;
for idel = 1:nDelta
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        subplot(noff_all,nDelta,((ioff-1)*nDelta)+idel)
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'ok')
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'tuning_byPref_byInt_sep_summary.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idel = 1:nDelta
    subplot(n,n2,idel)
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'o', 'Color', c(ioff,:));
        hold on
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
    xlabel('Stimulus Orientation (deg)')
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'tuning_byPref_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    for ioff = 1:noff_all
        for idel = 1:nDelta
            ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
            subplot(n,n2,idelta)
            errorbar(deltas(idel), mean(delta_resp_norm_all(ind_cells,idelta,ioff),1),std(delta_resp_norm_all(ind_cells,idelta,ioff),[],1)./sqrt(length(ind_cells)),'o', 'Color', c(ioff,:));
            hold on
        end
    end
    ylim([-0.1 1.2])
    xlabel('Cell preference group (deg)')
    title([num2str(deltas(idelta)) ' deg change'])
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'popTuning_byInt_summary.pdf'),'-dpdf','-fillpage')


%% Pref and OSI change
deltas = 22.5:22.5:180;
delta_diff = 180-deltas;
delta_diff(find(delta_diff>90)) = 180- delta_diff(find(delta_diff>90));
diffs = unique(delta_diff);
ndiff = length(diffs);

cell_group_n = zeros(1,ndiff);
group_list = zeros(size(max_dir_all,1),1);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    cell_group_n(1,idiff) = length(cell_ind);
    group_list(cell_ind,1) = idiff;
end

figure;
subplot(2,2,1)
osi_diff = zeros(noff, ndiff,2);
osi_avg = zeros(ndiff,2);
p_osi_diff = zeros(noff,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    osi_avg(idiff,1) = mean(OSI_all(cell_ind,3),1);
    osi_avg(idiff,2) = std(OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind));
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),1), std(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
        osi_diff(ioff,idiff,1) = mean(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),1);
        osi_diff(ioff,idiff,2) = std(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind));
        [h, p_osi_diff(ioff,idiff)] = ttest(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3));
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI'])
title('OSI')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,2)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_all(cell_ind,ioff)./OSI_all(cell_ind,3),1), std(OSI_all(cell_ind,ioff)./OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Ratio of OSI'])
title('OSI')
ylim([0.5 1.5])
xlim([-10 100])
hline(1)
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_k_all(cell_ind,ioff)-OSI_k_all(cell_ind,3),1), std(OSI_k_all(cell_ind,ioff)-OSI_k_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI-k'])
title('OSI-k')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,4)
HWHM_deg_all = rad2deg(HWHM_all);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(HWHM_deg_all(cell_ind,ioff)-HWHM_deg_all(cell_ind,3),1), std(HWHM_deg_all(cell_ind,ioff)-HWHM_deg_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in HWHM'])
title('HWHM')
ylim([-30 30])
xlim([-10 100])
hline(0)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- OSI by Interval'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'OSI_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
pref_ori_all_diff = pref_ori_all;
pref_ori_all_diff(find(pref_ori_all_diff>90)) = 180-pref_ori_all_diff(find(pref_ori_all_diff>90));
subplot(2,2,1)
g = jet(180);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all(iC, 2),pref_ori_all(iC, 1), 'o', 'Color', g(ceil(pref_ori_all(iC,3))+1,:))
    hold on
end
refline(1,0)
axis square
xlabel(['Pref Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Pref Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,2)
g = jet(91);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all_diff(iC, 2),pref_ori_all_diff(iC, 1), 'o', 'Color', g(ceil(pref_ori_all_diff(iC,3)+1),:))
    hold on
end
refline(1,0)
axis square
xlabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
pref_ori_diff = zeros(noff, ndiff,2);
p_pref_ori_diff = zeros(noff,ndiff);
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1), std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
        pref_ori_diff(ioff,idiff,1) = mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1);
        pref_ori_diff(ioff,idiff,2) = std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind));
        [h, p_pref_ori_diff(ioff,idiff)] = ttest(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3));
    end
end
xlabel(['Diff of Max from Adaptor (deg)'])
ylabel(['Change in Pref'])
ylim([-40 40])
hline(0)
xlim([-10 100])
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Change in Preference by Interval'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'prefDiff_byInt_summary.pdf'),'-dpdf','-fillpage')

col_mat = strvcat('b','r','y');
figure;
delta_resp_all_avg = zeros(noff_all, ndiff,2);
delta_resp_norm_all_avg = zeros(noff_all, ndiff,2);
p_delta_resp_norm = zeros(noff, ndiff);
for i = 1:noff_all
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        subplot(2,2,1)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        delta_resp_all_avg(i,idiff,1) = squeeze(mean(mean(delta_resp_all(good_ind_theta,del,i),2),1));
        delta_resp_all_avg(i,idiff,2) = squeeze(std(mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta));
        hold on
        subplot(2,2,2)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
        delta_resp_norm_all_avg(i,idiff,1) = squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1));
        delta_resp_norm_all_avg(i,idiff,2) = squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta));
        if i < 3
            [h p_delta_resp_norm(i,idiff)] = ttest(mean(delta_resp_norm_all(good_ind_theta,del,3),2), mean(delta_resp_norm_all(good_ind_theta,del,i),2));
        end
    end
end
subplot(2,2,1)
title('Absolute')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 .3])
subplot(2,2,2)
title('Normalized')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 1])
subplot(2,2,3)
delta_resp_group_avg = zeros(noff, ndiff,2);
p_delta_resp_group = zeros(noff,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    diff_resp_group= [];
    for i = 1:length(del)
        cell_ind = intersect(good_ind_theta, find(max_dir_all == del(i)));
        diff_resp_group = [diff_resp_group; squeeze(delta_resp_norm_all(cell_ind,del(i),:))];
    end
    for i = 1:noff
        errorbar(diffs(idiff), mean(diff_resp_group(:,i),1), std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),['o' col_mat(i,:)])
        hold on
        delta_resp_group_avg(i,idiff,1) = mean(diff_resp_group(:,i),1);
        delta_resp_group_avg(i,idiff,2) = std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1));
        [h p_delta_resp_group(i,idiff)] = ttest(diff_resp_group(:,3),diff_resp_group(:,i));
    end
end
title('Normalized')
ylabel('dF/F at max ori')
xlabel('Diff of Max from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Mean Resp by Interval'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'resp_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
subplot(2,2,1)
nCells = length(good_ind_theta);
delta_resp_sum = zeros(noff*nCells,ndiff);
start = 1;
off_i = [];
diff_i = [];
for i = 1:noff
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        delta_resp_sum(start:start+nCells-1,idiff) = mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2);
        hold on
    end
    start = start+nCells;
end
ylabel('Diff from control resp')
xlabel('Diff of Stim from Adaptor (deg)')
title('Abs')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 0.2])
subplot(2,2,2)
delta_resp_norm_sum = zeros(noff*nCells,ndiff);
start = 1;
off_ii = [];
diff_ii = [];
for i = 1:noff
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
        delta_resp_norm_sum(start:start+nCells-1,idiff) = mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2);
    end
    start = start+nCells;
end
ylabel('Diff from control resp')
xlabel('Diff of Stim from Adaptor (deg)')
title('Norm')
hline(0)
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([-0.2 0.4])

[n bin] = histc(pref_ori_all_diff(:,3),[0 20 70 91]);
figure
%subplot(2,2,3)
delta_resp_group = nan(sum(noff.*n,1),1);
start = 1;
off_id = [];
diff_id = [];
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    diff_resp_group= zeros(length(ind),noff);
    for i = 1:length(ind)
        [peak_val peak_loc] = max(fit_all(:,ind(i),3),[],1);
        fit_norm = squeeze(fit_all(:,ind(i),:))./peak_val;
        diff_resp_group(i,:) = fit_norm(peak_loc,1:2);
    end
    for i = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(diff_resp_group(:,i),1), std(pref_ori_all_diff(ind,3),[],1)./sqrt(size(diff_resp_group,1)),std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),{['o' col_mat(i,:)], col_mat(i,:),col_mat(i,:)})
        hold on
        delta_resp_group(start:start+length(ind)-1,1) = diff_resp_group(:,i);
        off_id = [off_id ones(1,length(ind))*offs(i)];
        diff_id = [diff_id ones(1,length(ind))*diffs(idiff)];
        start = start+length(ind);
    end
end
title('Norm')
ylabel('dF/F at max ori')
xlabel('Diff of Pref from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Diff Resp by Stimulus'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'diffResp_byStim_byFit_newBin2_summary.pdf'),'-dpdf','-fillpage')

[p_sum, table_sum, stats_sum] = anova2(delta_resp_sum, nCells);
[p_normsum, table_normsum, stats_normsum] = anova2(delta_resp_norm_sum, nCells);
[p_normgroup, table_normgroup, stats_normgroup] = anovan(delta_resp_group, {off_id, diff_id});

figure;
subplot(2,2,1)
pref_ori_group = nan(sum(noff.*n,1),1);
start = 1;
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    for ioff = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3),1), std(pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), std(pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), {['o' col_mat(ioff,:)], col_mat(ioff,:),col_mat(ioff,:)})
        hold on
        pref_ori_group(start:start+length(ind)-1,1) = pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3);
        start= start +length(ind);
    end
end
xlabel(['Diff of Pref from Adaptor (deg)'])
ylabel(['Difference in Pref (deg)'])
title('Pref')
ylim([-20 50])
hline(0)
xlim([-10 100])
subplot(2,2,2)
OSI_group = nan(sum(noff.*n,1),1);
start = 1;
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    for ioff = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(OSI_all(ind,ioff)-OSI_all(ind,3),1),std(pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), std(OSI_all(ind,ioff)-OSI_all(ind,3),[],1)./sqrt(length(ind)), {['o' col_mat(ioff,:)], col_mat(ioff,:),col_mat(ioff,:)})
        hold on
        OSI_group(start:start+length(ind)-1,1) = OSI_all(ind,ioff)-OSI_all(ind,3);
        start= start +length(ind);
    end
end
xlabel(['Diff of Pref from adaptor (deg)'])
ylabel(['Difference in OSI'])
title('OSI')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Change in Preference by Interval'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'diff_byInt_byPref_newBin_summary.pdf'),'-dpdf','-fillpage')
[p_prefdiff, table_prefdiff, stats_prefdiff] = anovan(pref_ori_group, {off_id, diff_id});
[p_osidiff, table_osidiff, stats_osidiff] = anovan(OSI_group, {off_id, diff_id});
%% max log likelihood

nboot = 1000;
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);
sz_fit = size(fit_all,1);

loglike_all_fun = nan(sz_fit,nDelta,noff_all,nboot+1);
loglike_all_fun_nosub = nan(sz_fit,nDelta,noff_all,nboot+1);
loglike_all_fun_subfact = nan(sz_fit,nDelta,noff_all,nboot+1);
maxloglike_all = nan(nDelta,noff_all,nboot+1);
maxloglike_all_subfact = nan(nDelta,noff_all,nboot+1);
maxloglike_all_nosub = nan(nDelta,noff_all,nboot+1);

[max_fit, max_ind] = max(fit_all(:,:,3),[],1);
delta_resp_all_thresh = bsxfun(@rdivide,delta_resp_all,max_fit').*100;
delta_resp_all_thresh(find(delta_resp_all_thresh<1)) = 1;
delta_resp_all_thresh(find(delta_resp_all_thresh>170)) = 170;
fit_all_thresh =  bsxfun(@rdivide,fit_all(:,:,3),max_fit)*100;
fit_all_thresh(find(fit_all_thresh<1)) = 1;
fit_all_thresh(find(fit_all_thresh>170)) = 170;
% [n_fit bin_fit] = histc(pref_ori_all(:,3), [0:22.5:180]);
% delta_avg = zeros(length(n_fit)-1, nDelta,noff_all);
% fit_avg = zeros(sz_fit, length(n_fit)-1);
% delta_avg_thresh = zeros(length(n_fit)-1, nDelta,noff_all);
% fit_avg_thresh = zeros(sz_fit, length(n_fit)-1);
% loglike_all_neurons = nan(length(n_fit)-1,sz_fit,nDelta,noff_all);
% loglike_all_fact = nan(length(n_fit)-1,1,nDelta,noff_all);
% loglike_all_sum = nan(length(n_fit)-1,sz_fit,nDelta,noff_all);
% for i = 1:length(n_fit)-1
%     ind = find(bin_fit == i);
%     fit_avg(:,i) = squeeze(mean(fit_all(:,ind,3),2))';
%     delta_avg(i,:,:) = mean(delta_resp_all(ind,:,:),1);
%     max_fit = max(fit_avg(:,i),[],1);
%     fit_avg_thresh(:,i) = (fit_avg(:,i)./max_fit)*100;
%     fit_avg_thresh(find(fit_avg_thresh<1)) = 1;
%     fit_avg_thresh(find(fit_avg_thresh>170)) = 170;
%     delta_avg_thresh(i,:,:) = (delta_avg(i,:,:)./max_fit)*100;
%     delta_avg_thresh(find(delta_avg_thresh<1)) = 1;
%     delta_avg_thresh(find(delta_avg_thresh>170)) = 170;
%     for idelta = 1:nDelta
%         for ioff =1:noff_all
%             loglike_all_neurons(i,:,idelta,ioff) = squeeze(log10(fit_avg_thresh(:,i)').*delta_avg_thresh(i,idelta,ioff))';
%             loglike_all_sum(i,:,idelta,ioff) = squeeze(fit_avg_thresh(:,i)');
%             loglike_all_fact(i,1,idelta,ioff) = log10(gamma(delta_avg_thresh(i,idelta,ioff)+1));
%         end
%     end
% end
% loglike_all_fun =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons,1), bsxfun(@plus,(nansum(loglike_all_sum,1)),nansum(loglike_all_fact,1))));
% loglike_all_fun_nosub =  squeeze(nansum(loglike_all_neurons,1));
% [max_val maxloglike_all] = max(loglike_all_fun,[],1);
% [max_val maxloglike_all_nosub] = max(loglike_all_fun_nosub,[],1);
% maxloglike_all = squeeze(maxloglike_all);
% maxloglike_all_nosub = squeeze(maxloglike_all_nosub);
loglike_all_neurons = nan(cellN,sz_fit,nDelta,noff_all);
loglike_all_fact = nan(cellN,1,nDelta,noff_all);
loglike_all_sum = nan(cellN,sz_fit,nDelta,noff_all);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    for idelta = 1:nDelta
        for ioff =1:noff_all
            loglike_all_neurons(iC,:,idelta,ioff) = squeeze(log10(fit_all_thresh(:,iC)').*delta_resp_all_thresh(iC,idelta,ioff))';
            loglike_all_sum(iC,:,idelta,ioff) = squeeze(fit_all_thresh(:,iC)');
            loglike_all_fact(iC,1,idelta,ioff) = log10(gamma(delta_resp_all_thresh(iC,idelta,ioff)+1));
        end
    end
end

%find best subtraction to account for inhomogeneity
loglike_all_fish = nan(sz_fit,nDelta,100);
maxloglike_fish = nan(nDelta,100);
for n = 1:100;
    loglike_all_fish(:,:,:,n) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(OSI_cells,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(OSI_cells,:,:,:),1)./n),nansum(loglike_all_fact(OSI_cells,:,:,:),1))));
    for idelta = 1:nDelta
        [max_val maxloglike_fish(idelta,n)] = max(loglike_all_fish(:,idelta,3,n),[],1);
    end
end
maxloglike_temp = maxloglike_fish;
ind = find(maxloglike_temp(8,:)<90);
maxloglike_temp(8,ind) = maxloglike_temp(8,ind)+180;
maxloglike_diff = bsxfun(@minus,maxloglike_temp, deltas');
maxloglike_sos = sum(maxloglike_diff.^2,1);
[min_val, min_ind] = min(maxloglike_sos,[],2);

for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = OSI_cells(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = OSI_cells;
    end
    loglike_all_fun(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(ind_cells_temp,:,:,:),1)./min_ind),nansum(loglike_all_fact(ind_cells_temp,:,:,:),1))));
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            [max_val maxloglike_all(idelta,ioff,iboot)] = max(loglike_all_fun(:,idelta,ioff,iboot),[],1);
        end
    end
end

% figure;
% for idelta = 1:nDelta
%     subplot(3,3,idelta)
%     for ioff = 1:noff_all
%         plot([0:180], loglike_all_fun(:,idelta,ioff));
%         hold on
%     end
%     title([num2str(deltas(idelta)) 'deg'])
%     xlabel('Orientation')
%     ylabel('Log likelihood')
%     xlim([-10 190])
%     %ylim([-3000 0])
% end

figure;
for idelta = 1:nDelta
    subplot(3,3,idelta)
    for ioff = 1:noff_all
        shadedErrorBar([0:180], loglike_all_fun(:,idelta,ioff,1), squeeze(std(loglike_all_fun(:,idelta,ioff,2:end),[],4)), {[col_mat(ioff) '-'],'markerfacecolor',col_mat(ioff)});
        hold on
    end
    title([num2str(deltas(idelta)) 'deg'])
    xlabel('Orientation')
    ylabel('Log likelihood')
    xlim([-10 190])
    %ylim([-3000 0])
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood Function- all fits- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_byDelta_byInt_summary_allFits.pdf'),'-dpdf','-fillpage')



%likelihood is 0
figure;
for ioff = 1:noff_all
    temp_likely = squeeze(loglike_all_fun(:,:,ioff,1));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun(end,:,ioff,2:end),[],4))]);
    hold on
    xlim([-22 202])
    xlabel('Orientation')
    ylabel('Log likelihood of 0 deg')
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood of 0 deg stimulus- all cells- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_is0_combo_allCells_100X.pdf'),'-dpdf','-fillpage')

figure;
for ioff = 1:noff_all
    subplot(2,1,1)
    plot([0 deltas], [maxloglike_all(end,ioff,1);maxloglike_all(:,ioff,1)],'-')
    hold on
    subplot(2,1,2)
    scatter(reshape(repmat([0 deltas], [nboot, 1])', [1, nboot*(nDelta+1)]), reshape(squeeze([maxloglike_all(end,ioff,2:nboot+1);maxloglike_all(:,ioff,2:nboot+1)]), [1, nboot*(nDelta+1)]),'o')
    hold on
end
for i = 1:2
    subplot(2,1,i)
    xlabel('Orientation')
    ylabel('Max log likelihood')
    axis square
end
subplot(2,1,1)
title('Orig')
subplot(2,1,2)
title('Bootstrap')
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Max Log Likelihood - scaled sub- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_scaledsub_100X.pdf'),'-dpdf','-fillpage')

maxloglike_change = maxloglike_all;
maxloglike_change(find(maxloglike_all>90)) = 179-maxloglike_change(find(maxloglike_all>90));
maxloglike_change_avg = zeros(ndiff,noff_all,2);
maxloglike_change_all = zeros(ndiff*nboot, noff);
start = 1;
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    maxloglike_change_avg(idiff,:,1) = mean(maxloglike_change(del,:,1),1);
    maxloglike_change_avg(idiff,:,2) = std(mean(maxloglike_change(del,:,2:nboot+1),1),[],3);
    maxloglike_change_all(start:start+nboot-1,:) = squeeze(mean(maxloglike_change(del,1:2,2:nboot+1),1))';
    start = start+nboot;
end
figure;
for i = 1:noff_all
	errorbar(diffs, maxloglike_change_avg(:,i,1), maxloglike_change_avg(:,i,2),'-o');
    hold on
end
xlabel('Stim - Adapter')
ylabel('Max likely ori - Adapter')
xlim([-10 100])
ylim([-10 100])
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Max Log Likelihood dist from Adapter - scaled sub- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_distAdapt_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_max, table_max, stats_max] = anova2(maxloglike_change_all, nboot);


figure;
delta_sq_boot = zeros(nDelta,nDelta,noff_all);
delta_sq = zeros(nDelta,nDelta,noff_all);
for ioff = 1:noff_all
    for idelta = 1:nDelta
        for idel = 1:nDelta
            delta_sq_boot(idel,idelta,ioff) = length(find(maxloglike_all(idelta,ioff,2:end) == idel))./1000;
        end
        delta_sq(maxloglike_all(idelta,ioff,1),idelta,ioff) = 1;
    end
    subplot(3,2,(ioff*2)-1)
    imagesc(flipud(delta_sq(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI'])
    colormap(hot)
    subplot(3,2,(ioff*2))
    imagesc(flipud(delta_sq_boot(:,:,ioff)),[0 1])
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI- bootstrap'])
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    colormap(hot)
    colorbar
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Max Log Likelihood- All cells- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_byInt_summary_allCells.pdf'),'-dpdf','-fillpage')

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxLogLike_summary.mat'),'delta_sq_boot','delta_sq','maxloglike_all','loglike_all_neurons','loglike_all_sum', 'loglike_all_fact')

%% max log likelihood by trials
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);
sz_fit = size(fit_all,1);
n_ind = zeros(nexp, nDelta,noff_all);
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            n_ind(iexp,idelta,ioff) = size(ppTemp{ioff, idelta},2);
        end
    end
end
min_trials = squeeze(min(n_ind,[],1));

[max_fit, max_ind] = max(fit_all(:,:,3),[],1);
fit_all_thresh =  bsxfun(@rdivide,fit_all(:,:,3),max_fit)*100;
fit_all_thresh(find(fit_all_thresh<1)) = 1;
fit_all_thresh(find(fit_all_thresh>170)) = 170;

ppResp_trials = cell(noff_all,nDelta);
start = 1;
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    nc = size(ppTemp{1,1},1);
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            ppT = ppTemp{ioff,idelta};
            ppT_norm = bsxfun(@rdivide, ppT(:,1:min_trials(idelta,ioff)), max_fit(1,start:nc+start-1)')*100;
            ppT_norm(find(ppT_norm<1)) = 1;
            ppT_norm(find(ppT_norm>170)) = 170;
            ppResp_trials{ioff,idelta}(start:nc+start-1,:) = ppT_norm;
        end
    end
    start = start+nc;
end

cellN = length(good_ind_theta);
nCells = size(fit_all_thresh,2);
nboot = 1000;
sz_fit = size(fit_all,1);
loglike_all_neurons = cell(nDelta,noff_all,nboot+1);
loglike_all_fact = cell(nDelta,noff_all,nboot+1);
loglike_all_sum = cell(nDelta,noff_all,nboot+1);
loglike_all_fun = cell(nDelta,noff_all,nboot+1);
maxloglike_all = cell(nDelta,noff_all,nboot+1);
maxloglike_change = cell(nDelta,noff_all,nboot+1);
maxloglike_change_diff = cell(ndiff,noff_all,nboot+1);
for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = good_ind_theta(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = good_ind_theta;
    end
    fprintf('.')
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            loglike_all_neurons{idelta,ioff,iboot} = zeros(sz_fit,cellN,min_trials(idelta,ioff));
            loglike_all_sum{idelta,ioff,iboot} = zeros(sz_fit,cellN,min_trials(idelta,ioff));
            loglike_all_fact{idelta,ioff,iboot} = zeros(1,cellN,min_trials(idelta,ioff));
            for i = 1:min_trials(idelta,ioff)
                for iCell = 1:cellN
                    iC = ind_cells_temp(iCell);
                    loglike_all_neurons{idelta,ioff,iboot}(:,iC,i) = squeeze(bsxfun(@times,log(fit_all_thresh(:,iC)'),ppResp_trials{ioff,idelta}(iC,i)))';
                    loglike_all_sum{idelta,ioff,iboot}(:,iC,i) = fit_all_thresh(:,iC);
                    loglike_all_fact{idelta,ioff,iboot}(1,iC,i) = log(gamma(ppResp_trials{ioff,idelta}(iC,i)+1));
                end
            end
            loglike_all_fun{idelta,ioff,iboot} =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons{idelta,ioff,iboot}(:,ind_cells_temp,:),2), bsxfun(@plus,(nansum(loglike_all_sum{idelta,ioff,iboot}(:,ind_cells_temp,:),2)),nansum(loglike_all_fact{idelta,ioff,iboot}(:,ind_cells_temp,:),2))));
            [max_val maxloglike_all{idelta,ioff,iboot}] = max(loglike_all_fun{idelta,ioff,iboot},[],1);
            maxloglike_change{idelta,ioff,iboot} = maxloglike_all{idelta,ioff,iboot};
            maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90)) = 179-maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90));
        end
    end

    if iboot == 0
        figure;
        for idelta = 1:nDelta
            for ioff = 1:noff_all
                subplot(3,3,idelta)
                shadedErrorBar(0:180, mean(loglike_all_fun{idelta,ioff,iboot},2), std(loglike_all_fun{idelta,ioff,iboot},[],2)./sqrt(size(loglike_all_fun{idelta,ioff,iboot},2)),{[col_mat(ioff) '-'],'markerfacecolor',col_mat(ioff)});
                hold on
            end
            title([num2str(deltas(idelta)) 'deg'])
            xlabel('Orientation')
            ylabel('Log likelihood')
            xlim([-10 190])
            %ylim([-3000 0])
        end
        suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood Function- all fits- by trial- 100X'])
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_byDelta_byInt_summary_allFits_byTrial.pdf'),'-dpdf','-fillpage')

        %likelihood is 0
        figure;
        likely_all = [];
        off_id = [];
        for ioff = 1:noff_all
            for idelta = 1:nDelta
                temp_likely = loglike_all_fun{idelta,ioff,iboot};
                errorbar(deltas(idelta), mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                hold on
                if idelta == nDelta
                    errorbar(0, mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                end
                if idelta == 1
                    likely_all = [likely_all temp_likely(end,:)];
                    off_id = [off_id ioff*ones(size(temp_likely(1,:)))];
                end
            end
            xlim([-22 202])
            xlabel('Orientation')
            ylabel('Log likelihood of 0 deg')
        end
        suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood of 0 deg stimulus- all cells- by trial- 100X'])
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_is0_combo_allCells_byTrial_100X.pdf'),'-dpdf','-fillpage')
        [p_like0, table_like0, stats_like0] = anovan(likely_all,{off_id});
    end

    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        for ioff = 1:noff_all
            for i =1:length(del)
                maxloglike_change_diff{idiff,ioff,iboot} = [maxloglike_change_diff{idiff,ioff,iboot} maxloglike_change{del(i),ioff,iboot}];
            end
        end
    end
end
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_change.mat'),'maxloglike_change_diff','maxloglike_change','loglike_all_fun')

figure;
maxloglike_change_diff_all = [];
maxloglike_change_diff_all_1 = [];
maxloglike_change_diff_all_2 = [];
off_id = [];
diff_id = [];
off_id_1 = [];
diff_id_1 = [];
off_id_2 = [];
diff_id_2 = [];
for idiff = 1:ndiff
    for ioff = 1:noff_all
        %errorbar(diffs(idiff),mean(maxloglike_change_diff{idiff,ioff,1},2), std(maxloglike_change_diff{idiff,ioff,1},[],2)./sqrt(size(maxloglike_change_diff{idiff,ioff,1},2)),['-o' col_mat(ioff,:)]);
        %hold on
        sz = size(maxloglike_change_diff{idiff,ioff,1});
        maxloglike_change_diff_all = [maxloglike_change_diff_all maxloglike_change_diff{idiff,ioff,1}];
        off_id = [off_id ioff*ones(sz)];
        diff_id = [diff_id idiff*ones(sz)];
        if idiff == 1
            off_id_1 = [off_id_1 ioff*ones(sz)];
            diff_id_1 = [diff_id_1 idiff*ones(sz)];
            maxloglike_change_diff_all_1 = [maxloglike_change_diff_all_1 maxloglike_change_diff{idiff,ioff,1}];
        end
        if idiff == 2
            off_id_2 = [off_id_2 ioff*ones(sz)];
            diff_id_2 = [diff_id_2 idiff*ones(sz)];
            maxloglike_change_diff_all_2 = [maxloglike_change_diff_all_2 maxloglike_change_diff{idiff,ioff,1}];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Max likely ori - Adapter')
xlim([-10 100])
ylim([-10 100])
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Max Log Likelihood dist from Adapter - by trial- scaled sub- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_max, table_max, stats_max] = anovan(maxloglike_change_diff_all,{off_id, diff_id});
[p_max_1, table_max_1, stats_max_1] = anovan(maxloglike_change_diff_all_1,{off_id_1});
[p_max_2, table_max_2, stats_max_2] = anovan(maxloglike_change_diff_all_2,{off_id_2});
[comp_max_1] = multcompare(stats_max_1);
[comp_max_2] = multcompare(stats_max_2);

figure;
loglike_diff_all = [];
loglike_diff_all_1 = [];
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    for ioff = 1:noff_all
        loglike_diff{idiff,ioff} = [];
        for idel = 1:length(del)
            loglike_diff{idiff,ioff} = [loglike_diff{idiff,ioff} loglike_all_fun{del(idel),ioff,1}];
        end
        temp_likely = loglike_diff{idiff,ioff};
        errorbar(diffs(idiff),mean(temp_likely(end,:),2), std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['-o' col_mat(ioff,:)]);
        hold on
        sz = size(loglike_diff{idiff,ioff});
        loglike_diff_all = [loglike_diff_all temp_likely(end,:)];
        if idiff == 1
            loglike_diff_all_1 = [loglike_diff_all_1 temp_likely(end,:)];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Likelyhood = 0')
xlim([-10 100])
ylim([-10000 0])
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood is 0 - by trial- scaled sub- 100X'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_is0_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_like0, table_like0, stats_like0] = anovan(loglike_diff_all,{off_id, diff_id});
[p_like0_1, table_like0_1, stats_like0_1] = anovan(loglike_diff_all_1,{off_id_1});
[comp_like0_1] = multcompare(stats_like0_1);
%save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_params.mat'),'loglike_all_neurons','loglike_all_sum','loglike_all_fact')

%% ROC
h = zeros(nDelta,1);
p = zeros(nDelta,1);
for idelta = 1:nDelta
    [h(idelta,:),p(idelta,:)] = ttest(roc_resp_all(:,1,idelta),roc_resp_all(:,2,idelta));
end
h_roc = zeros(nDelta,nDelta);
p_roc = zeros(nDelta,nDelta);
for i = 1:nDelta
    for idel = 1:nDelta
        ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
        [h_roc(idel,i), p_roc(idel,i)] = ttest(roc_resp_all(ind_cells,1,i),roc_resp_all(ind_cells,2,i));
    end
end

h_roc_diff = zeros(ndiff,ndiff);
p_roc_diff = zeros(ndiff,ndiff);
for i = 1:ndiff
    idiff = find(delta_diff == diffs(i)); 
    for idel = 1:ndiff
        dels = find(delta_diff == diffs(idel));
        ind_cells = [];
        for ii = 1:length(dels)
            ind_cells = [ind_cells; intersect(good_ind_theta, find(max_dir_all == dels(ii)))];
        end
        [h_roc_diff(idel,i), p_roc_diff(idel,i)] = ttest(mean(roc_resp_all(ind_cells,1,idiff),3),mean(roc_resp_all(ind_cells,2,idiff),3));
    end
end

figure
avg_roc = zeros(2,2,ndiff);
roc_all = [];
idiff_all = [];
ioff_all = [];
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(mean(roc_resp_all(good_ind_theta,1,del),3), mean(roc_resp_all(good_ind_theta,2,del),3), 'ob')
    roc_all = [roc_all ([mean(roc_resp_all(good_ind_theta,1,del),3); mean(roc_resp_all(good_ind_theta,2,del),3)])];
    avg_roc(:,1,idiff) = [mean(mean(roc_resp_all(good_ind_theta,1,del),3),1); mean(mean(roc_resp_all(good_ind_theta,2,del),3),1)];
    avg_roc(:,2,idiff) = [std(mean(roc_resp_all(good_ind_theta,1,del),3),[],1)./sqrt(length(good_ind_theta)); std(mean(roc_resp_all(good_ind_theta,2,del),3),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(mean(roc_resp_all(good_ind_theta,1,del),3), mean(roc_resp_all(good_ind_theta,2,del),3));
    axis square
    xlim([0 1])
    ylim([0 1])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc(:,1,:))', squeeze(avg_roc(:,2,:))', '-o')
ylim([0.3 0.7])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- ROC 0 (short ISI) vs Diff from Adapter']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'rocv180_allCells_summary.pdf'),'-dpdf','-fillpage')
[p_roc_all, table_roc_all, stats_roc_all] = anova2(roc_all, length(good_ind_theta));
comp_roc_all = multcompare(stats_roc_all);

figure
avg_roc_abs = zeros(2,2,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5), abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5), 'ob')
    avg_roc_abs(:,1,idiff) = [mean(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5),1); mean(abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5),1)];
    avg_roc_abs(:,2,idiff) = [std(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5),[],1)./sqrt(length(good_ind_theta)); std(abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5), abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5));
    axis square
    xlim([0 .5])
    ylim([0 .5])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc_abs(:,1,:))', squeeze(avg_roc_abs(:,2,:))', '-o')
ylim([0 0.5])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- ROC 0 (short ISI) vs Diff from Adapter- Absolute value']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'rocv180_allCells_summary_abs.pdf'),'-dpdf','-fillpage')

%% tuning figures- sub
delta_resp_norm_all = bsxfun(@rdivide, delta_resp_sub_all, max(delta_resp_sub_all(:,:,3),[],2));
figure;
for idel = 1:nDelta
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        subplot(noff_all,nDelta,((ioff-1)*nDelta)+idel)
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'ok')
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized- sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'tuning_byPref_byInt_sep_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idel = 1:nDelta
    subplot(n,n2,idel)
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'o', 'Color', c(ioff,:));
        hold on
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
    xlabel('Stimulus Orientation (deg)')
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized- sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'tuning_byPref_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    for ioff = 1:noff_all
        for idel = 1:nDelta
            ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
            subplot(n,n2,idelta)
            errorbar(deltas(idel), mean(delta_resp_norm_all(ind_cells,idelta,ioff),1),std(delta_resp_norm_all(ind_cells,idelta,ioff),[],1)./sqrt(length(ind_cells)),'o', 'Color', c(ioff,:));
            hold on
        end
    end
    ylim([-0.1 1.2])
    xlabel('Cell preference group (deg)')
    title([num2str(deltas(idelta)) ' deg change'])
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- avg all cells- normalized- sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'popTuning_byInt_summary_sub.pdf'),'-dpdf','-fillpage')


%% Pref and OSI change - sub
deltas = 22.5:22.5:180;
delta_diff = 180-deltas;
delta_diff(find(delta_diff>90)) = 180- delta_diff(find(delta_diff>90));
diffs = unique(delta_diff);
ndiff = length(diffs);

figure;
subplot(2,2,1)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_sub_all(cell_ind,ioff)-OSI_sub_all(cell_ind,3),1), std(OSI_sub_all(cell_ind,ioff)-OSI_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI'])
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,2)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_sub_all(cell_ind,ioff)./OSI_sub_all(cell_ind,3),1), std(OSI_sub_all(cell_ind,ioff)./OSI_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Ratio of OSI'])
ylim([0.5 1.5])
xlim([-10 100])
hline(1)
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_k_sub_all(cell_ind,ioff)-OSI_k_sub_all(cell_ind,3),1), std(OSI_k_sub_all(cell_ind,ioff)-OSI_k_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI-k'])
title('OSI-k')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,4)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(HWHM_sub_all(cell_ind,ioff)-HWHM_sub_all(cell_ind,3),1), std(HWHM_sub_all(cell_ind,ioff)-HWHM_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in HWHM'])
title('HWHM')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- OSI by Interval- Sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'OSI_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
pref_ori_all_diff = pref_ori_sub_all;
pref_ori_all_diff(find(pref_ori_all_diff>90)) = 180-pref_ori_all_diff(find(pref_ori_all_diff>90));
subplot(2,2,1)
g = jet(180);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_sub_all(iC, 2),pref_ori_sub_all(iC, 1), 'o', 'Color', g(ceil(pref_ori_sub_all(iC,3))+1,:))
    hold on
end
refline(1,0)
axis square
xlabel(['Pref Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Pref Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,2)
g = jet(91);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all_diff(iC, 2),pref_ori_all_diff(iC, 1), 'o', 'Color', g(ceil(pref_ori_all_diff(iC,3)+1),:))
    hold on
end
refline(1,0)
axis square
xlabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1), std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of Max from Adaptor (deg)'])
ylabel(['Change in Pref'])
ylim([-40 40])
hline(0)
xlim([-10 100])
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Change in Preference by Interval- Sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'prefDiff_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

col_mat = strvcat('b','r','y');
figure;
for i = 1:noff_all
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        subplot(2,2,1)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_sub_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_sub_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
        subplot(2,2,2)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
    end
end
subplot(2,2,1)
title('Absolute')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 .3])
subplot(2,2,2)
title('Normalized')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 1])
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    diff_resp_group= [];
    for i = 1:length(del)
        cell_ind = intersect(good_ind_theta, find(max_dir_all == del(i)));
        diff_resp_group = [diff_resp_group; squeeze(delta_resp_norm_all(cell_ind,del(i),1:2))];
    end
    for i = 1:noff
        errorbar(diffs(idiff), mean(diff_resp_group(:,i),1), std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),['o' col_mat(i,:)])
        hold on
    end
end
title('Normalized')
ylabel('dF/F at max ori')
xlabel('Diff of Max from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Mean Resp by Interval- sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'resp_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

%% max log likelihood- sub

nboot = 1000;
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);

loglike_all_neurons = zeros(cellN,nDelta,nDelta,noff_all);
loglike_all_fact = zeros(cellN,1,nDelta,noff_all);
loglike_all_sum = zeros(cellN,nDelta,nDelta,noff_all);
loglike_all_fun = zeros(nDelta,nDelta,noff_all,nboot+1);
loglike_all_fun_nosub = zeros(nDelta,nDelta,noff_all,nboot+1);
loglike_all_fun_subfact = zeros(nDelta,nDelta,noff_all,nboot+1);
maxloglike_all = zeros(nDelta,noff_all,nboot+1);
maxloglike_all_subfact = zeros(nDelta,noff_all,nboot+1);
maxloglike_all_nosub = zeros(nDelta,noff_all,nboot+1);

delta_resp_all_thresh = delta_resp_sub_all*100;
delta_resp_all_thresh(find(delta_resp_sub_all<0)) = 0.001;
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    for idelta = 1:nDelta
        for ioff =1:noff_all
            loglike_all_neurons(iC,:,idelta,ioff) = squeeze(log10(delta_resp_all_thresh(iC,:,noff_all)).*delta_resp_all_thresh(iC,idelta,ioff))';
            loglike_all_sum(iC,:,idelta,ioff) = squeeze(delta_resp_all_thresh(iC,:,noff_all))';
            loglike_all_fact(iC,1,idelta,ioff) = log10(gamma(delta_resp_all_thresh(iC,idelta,ioff)+1));
        end
    end
end
    
for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = OSI_cells(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = OSI_cells;
    end

    loglike_all_fun(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(ind_cells_temp,:,:,:),1)),nansum(loglike_all_fact(ind_cells_temp,:,:,:),1))));
    loglike_all_fun_subfact(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), nansum(loglike_all_fact(ind_cells_temp,:,:,:),1)));
    loglike_all_fun_nosub(:,:,:,iboot) =  squeeze(nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1));
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            [max_val maxloglike_all(idelta,ioff,iboot)] = max(loglike_all_fun(:,idelta,ioff,iboot),[],1);
            [max_val maxloglike_all_subfact(idelta,ioff,iboot)] = max(loglike_all_fun_subfact(:,idelta,ioff,iboot),[],1);
            [max_val maxloglike_all_nosub(idelta,ioff,iboot)] = max(loglike_all_fun_nosub(:,idelta,ioff,iboot),[],1);
        end
    end
end


figure;
for idelta = 1:nDelta
    subplot(3,3,idelta)
    for ioff = 1:noff_all
        plot(deltas, loglike_all_fun(:,idelta,ioff,1))
        hold on
    end
    title([num2str(deltas(idelta)) 'deg'])
    xlabel('Orientation')
    ylabel('Log likelihood')
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood Function- all cells- 100X- Sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_byDelta_byInt_summary_allCells_100X_sub.pdf'),'-dpdf','-fillpage')



figure;
for ioff = 1:noff_all
    subplot(1,3,1)
    temp_likely = squeeze(loglike_all_fun_nosub(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun_nosub(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun_nosub(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun_nosub(end,:,ioff,2:end),[],4))]);
    hold on
    if ioff == noff_all
        title('No sub')
    end
    subplot(1,3,2)
    temp_likely = squeeze(loglike_all_fun(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun(end,:,ioff,2:end),[],4))]);
    
    hold on
    if ioff == noff_all
        title('Sub Fact and Sum')
    end
    
    subplot(1,3,3)
    temp_likely = squeeze(loglike_all_fun_subfact(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun_subfact(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun_subfact(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun_subfact(end,:,ioff,2:end),[],4))]);
    
    hold on
    if ioff == noff_all
        title('Sub Fact only')
    end
end
for i = 1:3
    subplot(1,3,i)
    xlim([-22 202])
    xlabel('Orientation')
    ylabel('Log likelihood of 0 deg')
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Log Likelihood of 0 deg stimulus- all cells- 100X- Sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'loglike_is0_combo_allCells_100X_sub.pdf'),'-dpdf','-fillpage')


figure;
delta_sq_boot = zeros(nDelta,nDelta,noff_all);
delta_sq = zeros(nDelta,nDelta,noff_all);
for ioff = 1:noff_all
    for idelta = 1:nDelta
        for idel = 1:nDelta
            delta_sq_boot(idel,idelta,ioff) = length(find(maxloglike_all(idelta,ioff,2:end) == idel))./1000;
        end
        delta_sq(maxloglike_all(idelta,ioff,1),idelta,ioff) = 1;
    end
    subplot(3,2,(ioff*2)-1)
    imagesc(flipud(delta_sq(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI'])
    colormap(hot)
    subplot(3,2,(ioff*2))
    imagesc(flipud(delta_sq_boot(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI- bootstrap'])
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    colormap(hot)
    colorbar
end
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- Max Log Likelihood- All cells- 100X- Sub'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'maxloglike_byInt_summary_allCells_100X_sub.pdf'),'-dpdf','-fillpage')

%% ROC - sub
h = zeros(nDelta,1);
p = zeros(nDelta,1);
for idelta = 1:nDelta
    [h(idelta,:),p(idelta,:)] = ttest(roc_resp_sub_all(:,1,idelta),roc_resp_sub_all(:,2,idelta));
end
h_roc = zeros(nDelta,nDelta);
p_roc = zeros(nDelta,nDelta);
for i = 1:nDelta
    for idel = 1:nDelta
        ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
        [h_roc(idel,i), p_roc(idel,i)] = ttest(roc_resp_sub_all(ind_cells,1,i),roc_resp_sub_all(ind_cells,2,i));
    end
end

figure
avg_roc = zeros(2,2,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(mean(roc_resp_sub_all(good_ind_theta,1,del),3), mean(roc_resp_sub_all(good_ind_theta,2,del),3), 'ob')
    avg_roc(:,1,idiff) = [mean(mean(roc_resp_sub_all(good_ind_theta,1,del),3),1); mean(mean(roc_resp_sub_all(good_ind_theta,2,del),3),1)];
    avg_roc(:,2,idiff) = [std(mean(roc_resp_sub_all(good_ind_theta,1,del),3),[],1)./sqrt(length(good_ind_theta)); std(mean(roc_resp_sub_all(good_ind_theta,2,del),3),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(mean(roc_resp_sub_all(good_ind_theta,1,del),3), mean(roc_resp_sub_all(good_ind_theta,2,del),3));
    axis square
    xlim([0 1])
    ylim([0 1])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc(:,1,:))', squeeze(avg_roc(:,2,:))', '-o')
ylim([0.3 0.7])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([reshape(flipud(rot90(mouse_mat)),[1 nexp*4]) '- ROC 0 (short ISI) vs Diff from Adapter - Sub']) 
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'rocv180_allCells_summary_sub.pdf'),'-dpdf','-fillpage')



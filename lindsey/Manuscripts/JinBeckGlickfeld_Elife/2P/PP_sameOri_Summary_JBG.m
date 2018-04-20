% current datasets
%% V1 all
dataset = 'V1all_';
mouse_mat = strvcat('i684', 'i689','i720','i738','i739');
date_mat = strvcat('170316', '170327','170802','170710','170802');
run_str_mat = strvcat('runs-002-004','runs-002-003','runs-002-003','runs-002-003','runs-002-003');

%%
nexp = size(mouse_mat,1);
mouse_str = [];
for iexp = 1:nexp
    mouse_str = [mouse_str mouse_mat(iexp,:)];
end

%%
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%%
ind_min = cell(1,nexp);
ind_min_alldir = cell(1,nexp);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat(iexp,:);
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    
    ind{iexp} = zeros(ndir,noff);
    for idir = 1:ndir
        for ioff = 1:noff
            ind{iexp}(idir,ioff) = length(intersect(find(tFramesOff(:,end) == offs(ioff)), find(baseDir == dirs(idir))));
        end
    end
    [ind_min_alldir{iexp} i] = min(sum(ind{iexp},1),[],2);
    [ind_min{iexp} i] = min(min(ind{iexp},[],1),[],2);
end

%% collect datasets

max_dir_all = [];
resp_all = [];
resp_all_sub = [];
resp_dir_all = [];
resp_dir_all_sub = [];
good_ind_temp = [];
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat(iexp,:);
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
    max_dir_all = [max_dir_all; max_dir];
    resp_all = cat(1,resp_all, resp_off);
    resp_all_sub = cat(1,resp_all_sub, resp_off_sub);
    resp_dir_all = cat(1,resp_dir_all, resp_off_dir);
    resp_dir_all_sub = cat(1,resp_dir_all_sub, resp_off_dir_sub);
    good_ind_ind = zeros(size(resp_off,1),1);
    good_ind_ind(good_ind) = 1;
    good_ind_temp = [good_ind_temp; good_ind_ind];
end
load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

good_ind_all = find(good_ind_temp);
resp_maxdir_all = nan(size(resp_dir_all,1), size(resp_dir_all,3));
resp_nextdir_all = nan(size(resp_dir_all,1), size(resp_dir_all,3));
resp_maxdir_all_sub = nan(size(resp_dir_all,1), size(resp_dir_all,3));
ndir = length(unique(max_dir_all));
for iCell = 1:length(good_ind_all)
    iC = good_ind_all(iCell);
    x =max_dir_all(iC,1);
    resp_maxdir_all(iC,:) = squeeze(resp_dir_all(iC,x,:));
    resp_maxdir_all_sub(iC,:) = squeeze(resp_dir_all_sub(iC,max_dir_all(iC,1),:));
    y = x+1;
    z = x-1;
    if y > ndir
        y = 1;
    end
    if z < 1
        z = ndir;
    end 
    resp_nextdir_all(iC,:) = squeeze(mean(resp_dir_all(iC,[y z],:),2));
end

norm_all = bsxfun(@rdivide, resp_all, resp_all(:,end));
norm_all_sub = bsxfun(@rdivide, resp_all_sub, resp_all_sub(:,end));
norm_nextdir_all = bsxfun(@rdivide, resp_nextdir_all, resp_nextdir_all(:,end));
norm_maxdir_all = bsxfun(@rdivide, resp_maxdir_all, resp_maxdir_all(:,end));
norm_maxdir_all_sub = bsxfun(@rdivide, resp_maxdir_all_sub, resp_maxdir_all_sub(:,end));
    
%% figures
%good_ind_all = intersect(good_ind_all, find(resp_all(:,6)>0.02));
norm_all_avg = mean(norm_all(good_ind_all,:),1);
norm_all_sem = std(norm_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all));
[norm_all_h, norm_all_p] = ttest(norm_all,1);
[p_norm_all table_norm_all stats_norm_all] = anova1(norm_all);
[comp_norm_all] = multcompare(stats_norm_all);
resp_ratio_fit = fitlm(resp_all(good_ind_all,6), norm_all(good_ind_all,1),'linear');

norm_maxdir_all_avg = mean(norm_maxdir_all(good_ind_all,:),1);
norm_maxdir_all_sem = std(norm_maxdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all));
norm_nextdir_all_avg = mean(norm_nextdir_all(good_ind_all,:),1);
norm_nextdir_all_sem = std(norm_nextdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all));
resp_maxdir_all_avg = mean(resp_maxdir_all(good_ind_all,:),1);
resp_maxdir_all_sem = std(resp_maxdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all));
resp_nextdir_all_avg = mean(resp_nextdir_all(good_ind_all,:),1);
resp_nextdir_all_sem = std(resp_nextdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all));

[norm_maxdir_all_h, norm_maxdir_all_p] = ttest(norm_maxdir_all,1);
[p_norm_maxdir_all table_norm_maxdir_all stats_norm_maxdir_all] = anova1(norm_maxdir_all);
[comp_norm_maxdir_all] = multcompare(stats_norm_maxdir_all);
resp_maxdir_ratio_fit = fitlm(resp_maxdir_all(good_ind_all,6), norm_maxdir_all(good_ind_all,1),'linear');

[norm_maxVnext_h, norm_maxVnext_p] = ttest(norm_maxdir_all(:,1), norm_nextdir_all(:,1));
[resp_maxVnext_h, resp_maxVnext_p] = ttest(resp_maxdir_all(:,6), resp_nextdir_all(:,6));

edges = [0 0.05 0.15 0.3 0.8];
[n,bin] = histc(resp_all(:,6),edges);
figure;
subplot(2,2,1)
for i = 1:length(n)-1
    avg_resp = mean(resp_all(intersect(good_ind_all,find(bin==i)),6),1);
    sem_resp = std(resp_all(intersect(good_ind_all,find(bin==i)),6),[],1)./sqrt(length(intersect(good_ind_all,find(bin==i))));
    avg_norm = mean(norm_all(intersect(good_ind_all,find(bin==i)),1),1);
    sem_norm = std(norm_all(intersect(good_ind_all,find(bin==i)),1),[],1)./sqrt(length(intersect(good_ind_all,find(bin==i))));
    errorbarxy(avg_resp,avg_norm,sem_resp,sem_norm);
    hold on
end
xlim([0 0.5])
xlabel('Response amplitude (dF/F)')
ylim([0 1])
ylabel('Normalized amplitude')
subplot(2,2,2)
errorbarxy(resp_maxdir_all_avg(:,6),norm_maxdir_all_avg(:,1),resp_maxdir_all_sem(:,6),norm_maxdir_all_sem(:,1))
hold on
errorbarxy(resp_nextdir_all_avg(:,6),norm_nextdir_all_avg(:,1),resp_nextdir_all_sem(:,6),norm_nextdir_all_sem(:,1))
xlim([0 0.2])
xlabel('Response amplitude (dF/F)')
ylim([0 1])
ylabel('Normalized amplitude')
legend('Pref Ori','Neighbor Ori')
suptitle(['Paired Pulse- Same Ori- ' dataset(1,1:length(dataset)-1) ' ' mouse_str ' - n = ' num2str(length(good_ind_all)) ' Cells'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation', [dataset 'ppAdaptation_controls.pdf']),'-dpdf','-fillpage')
[p] = anova1(norm_all(good_ind_all,1),bin(good_ind_all));

x_range = [0:240*(1000/frameRateHz)];
off_all= [offs; 240];
figure;
subplot(2,2,1)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_all(good_ind_all,:),1)),squeeze(std(norm_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_all(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');
title(['All dirs- No sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
subplot(2,2,2)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_maxdir_all(good_ind_all,:),1)),squeeze(std(norm_maxdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_maxdir_all(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');
title(['Max dirs only- No sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
subplot(2,2,3)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_all_sub(good_ind_all,:),1)),squeeze(std(norm_all_sub(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_all_sub(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');
title(['All dirs- Sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
subplot(2,2,4)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_maxdir_all_sub(good_ind_all,:),1)),squeeze(std(norm_maxdir_all_sub(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_maxdir_all_sub(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');title(['Max dirs only- Sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
suptitle(['Paired Pulse- Same Ori- ' dataset(1,1:length(dataset)-1) ' ' mouse_str ' - n = ' num2str(length(good_ind_all)) ' Cells'])

print(fullfile(LG_base, '\Analysis\2P', 'Adaptation', [dataset 'ppAdaptation_summary.pdf']),'-dpdf','-fillpage')
save(fullfile(LG_base, '\Analysis\2P', 'Adaptation', [dataset 'ppAdaptation_summary.mat']),'resp_all', 'norm_all_avg', 'norm_all_sem', 'mouse_mat', 'date_mat', 'run_str_mat', 'norm_all', 'norm_all_sub', 'norm_maxdir_all', 'norm_maxdir_all_sub', 'good_ind_all');

%% compare electrophysiology to imaging data

load(fullfile(LG_base, '\Analysis\2P', 'Adaptation', 'ppMiao.mat'))

[val, ind] = sort(resp_all(good_ind_all,6),'descend');
for i = 1:length(good_ind_all)
pp_all = [norm_all(good_ind_all(ind(1:i)),1:5); pp(:,1:5)];
pp_all_reshape = reshape(pp_all,[size(pp_all,1)*size(pp_all,2),1]);
g_all = repmat([ones(length(good_ind_all(ind(1:i))),1); 2.*ones(size(pp,1),1)],[5 1]);
stim_all = [ones(size(pp_all,1),1); 2.*ones(size(pp_all,1),1); 3.*ones(size(pp_all,1),1); 4.*ones(size(pp_all,1),1); 5.*ones(size(pp_all,1),1)];
p(i,:) = anovan(pp_all_reshape,{g_all,stim_all});
end

figure; plot(val)
hold on; plot(p(:,1))
hline(0.05)
vline(30)
legend('Amplitude','P')

pp_all = [norm_all(good_ind_all,1:5); pp(:,1:5)];
pp_all_reshape = reshape(pp_all,[size(pp_all,1)*size(pp_all,2),1]);
g_all = repmat([ones(length(good_ind_all),1); 2.*ones(size(pp,1),1)],[5 1]);
stim_all = [ones(size(pp_all,1),1); 2.*ones(size(pp_all,1),1); 3.*ones(size(pp_all,1),1); 4.*ones(size(pp_all,1),1); 5.*ones(size(pp_all,1),1)];
[p_IvsE, table_IvsE, stats_IvsE] = anovan(pp_all_reshape,{g_all,stim_all});
out = multcompare(stats_IvsE);



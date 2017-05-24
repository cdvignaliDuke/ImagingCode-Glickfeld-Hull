%% current datasets

mouse_mat = strvcat('i684', 'i689');
date_mat = strvcat('170316', '170327');
run_str_mat = strvcat('runs-002-004','runs-002-003');

nexp = size(mouse_mat,1);
mouse_str = [];
for iexp = 1:nexp
    mouse_str = [mouse_str mouse_mat(iexp,:)];
end

%%
ind_min = cell(1,nexp);
ind_min_alldir = cell(1,nexp);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat(iexp,:);
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    
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
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
    max_dir_all = [max_dir_all; max_dir];
    resp_all = cat(1,resp_all, resp_off);
    resp_all_sub = cat(1,resp_all_sub, resp_off_sub);
    resp_dir_all = cat(1,resp_dir_all, resp_off_dir);
    resp_dir_all_sub = cat(1,resp_dir_all_sub, resp_off_dir_sub);
    good_ind_ind = zeros(size(resp_off,1),1);
    good_ind_ind(good_ind) = 1;
    good_ind_temp = [good_ind_temp; good_ind_ind];
end
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

good_ind_all = find(good_ind_temp);
resp_maxdir_all = nan(size(resp_dir_all,1), size(resp_dir_all,3));
resp_maxdir_all_sub = nan(size(resp_dir_all,1), size(resp_dir_all,3));
for iCell = 1:length(good_ind_all)
    iC = good_ind_all(iCell);
    resp_maxdir_all(iC,:) = squeeze(resp_dir_all(iC,max_dir_all(iC,1),:));
    resp_maxdir_all_sub(iC,:) = squeeze(resp_dir_all_sub(iC,max_dir_all(iC,1),:));
end

norm_all = bsxfun(@rdivide, resp_all, resp_all(:,end));
norm_all_sub = bsxfun(@rdivide, resp_all_sub, resp_all_sub(:,end));
norm_maxdir_all = bsxfun(@rdivide, resp_maxdir_all, resp_maxdir_all(:,end));
norm_maxdir_all_sub = bsxfun(@rdivide, resp_maxdir_all_sub, resp_maxdir_all_sub(:,end));
    
%% figures
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
subplot(2,2,2)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_maxdir_all(good_ind_all,:),1)),squeeze(std(norm_maxdir_all(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_maxdir_all(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');
title(['Max dirs only- No sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
ylim([0 1.2])
subplot(2,2,3)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_all_sub(good_ind_all,:),1)),squeeze(std(norm_all_sub(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_all_sub(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');
title(['All dirs- Sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
subplot(2,2,4)
errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_maxdir_all_sub(good_ind_all,:),1)),squeeze(std(norm_maxdir_all_sub(good_ind_all,:),[],1)./sqrt(length(good_ind_all))), 'ok');
[tau, A,sse,R_square] = SingleExponFit(off_all.*(1000/frameRateHz), mean(norm_maxdir_all_sub(good_ind_all,:),1)');
fit_plot = 1-A.*exp(x_range./(-tau));
hold on; plot(x_range,fit_plot,'r');title(['Max dirs only- Sub- Tau-' num2str(chop(tau,3)) '; R^2- ' num2str(chop(R_square,2))])
ylim([0 1.2])
suptitle(['Paired Pulse- Same Ori- ' mouse_str ' - n = ' num2str(length(good_ind_all)) ' Cells'])

print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'ppAdaptation_summary.pdf'),'-dpdf','-fillpage')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'ppAdaptation_summary.mat'),'mouse_mat', 'date_mat', 'run_str_mat', 'norm_all', 'norm_all_sub', 'norm_maxdir_all', 'norm_maxdir_all_sub', 'good_ind_all');
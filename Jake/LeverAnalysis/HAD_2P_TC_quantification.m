% finds and maps responsive cells for each condition
% 1. calculate average timecourses for press/release events
% 2. calculate variability by trial over base and resp windows
% 3. ttest for significant responses
% 4. define cells by response to event

load([dest '_ROI_TCs.mat']);
load([dest_sub '_release_movies.mat'])
load([dest_sub '_press_movies.mat'])
nCells = size(success_movie,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n)
    n2= n-1;
else
    n2 = n;
end

%% 1. calculate average timecourses for press/release events
%average and sem across all release
release_movie = cat(1,success_movie, fail_movie);
avg_release = squeeze(mean(release_movie,1));
sem_release = squeeze(std(release_movie,1)./sqrt(size(release_movie,1)));
%average and sem across cells (and trials)
avg_release_all = mean(avg_release,1);
sem_release_all = std(avg_release,1)./sqrt(size(avg_release,1));

%average and sem across success
avg_success = squeeze(mean(success_movie,1));
sem_success = squeeze(std(success_movie,1)./sqrt(size(success_movie,1)));
%average and sem across cells (and trials)
avg_success_all = mean(avg_success,1);
sem_success_all = std(avg_success,1)./sqrt(size(avg_success,1));

%average and sem across failures
avg_fail = squeeze(mean(fail_movie,1));
sem_fail = squeeze(std(fail_movie,1)./sqrt(size(fail_movie,1)));
%average and sem across ROIs (and trials)
avg_fail_all = mean(avg_fail,1);
sem_fail_all = std(avg_fail,1)./sqrt(size(avg_fail,1));

%average and sem across presses
avg_press = squeeze(mean(press_long_movie,1));
sem_press = squeeze(std(press_long_movie,1)./sqrt(size(press_long_movie,1)));
%average and sem across ROIs (and trials)
avg_press_all = mean(avg_press,1);
sem_press_all = std(avg_press,1)./sqrt(size(avg_press,1));

save([dest_sub '_cell_TCs.mat'], 'avg_release', 'avg_press','avg_success', 'avg_fail', 'sem_release', 'sem_press','sem_success', 'sem_fail')
%% 2. calculate response amplitude and variablity by trial over base and resp windows
%by release
base_release_window = 1:pre_release_frames;
resp_release_window = pre_release_frames+1:pre_release_frames+round(100./double(ifi));

%release
release_base = squeeze(mean(release_movie(:,:,base_release_window),3));
release_resp = squeeze(mean(release_movie(:,:,resp_release_window),3));
release_resp_avg = mean((release_resp-release_base),1);
release_resp_sem = std((release_resp-release_base),[],1)./sqrt(size(release_resp,1));

%success
success_base = squeeze(mean(success_movie(:,:,base_release_window),3));
success_resp = squeeze(mean(success_movie(:,:,resp_release_window),3));
success_resp_avg = mean((success_resp-success_base),1);
success_resp_sem = std((success_resp-success_base),[],1)./sqrt(size(success_resp,1));

%fail
fail_base = squeeze(mean(fail_movie(:,:,base_release_window),3));
fail_resp = squeeze(mean(fail_movie(:,:,resp_release_window),3));
fail_resp_avg = mean((fail_resp-fail_base),1);
fail_resp_sem = std((fail_resp-fail_base),[],1)./sqrt(size(fail_resp,1));

%by press
base_press_window = 1:pre_press_frames-round(300./double(ifi));
resp_press_window = pre_press_frames-round(100./double(ifi)):pre_press_frames+round(100./double(ifi));
[maxx, indp] = max(mean(press_long_movie(:,:,resp_press_window),1),[],3);
press_base = zeros(size(press_long_movie,1),nCells);
press_resp = zeros(size(press_long_movie,1),nCells);
for ic = 1:nCells
    press_base(:,ic) = squeeze(mean(press_long_movie(:,ic,base_press_window),3));
    press_resp(:,ic) = squeeze(mean(press_long_movie(:,ic,resp_press_window(indp(ic))-round(50./double(ifi)):resp_press_window(indp(ic))+round(50./double(ifi))),3));
end
press_resp_avg = mean((press_resp-press_base),1);
press_resp_sem = std((press_resp-press_base),[],1)./sqrt(size(press_resp,1));

save([dest_sub '_cell_resp.mat'], 'base_release_window', 'resp_release_window', 'base_press_window', 'resp_press_window', 'press_base', 'press_resp', 'success_base', 'success_resp', 'fail_base', 'fail_resp', 'release_resp', 'release_base');

%% 2. ttest for significant responses
[release_h, release_p] = ttest(release_base, release_resp, 'dim', 1, 'tail', 'both');
[success_h, success_p] = ttest(success_base, success_resp, 'dim', 1, 'tail', 'both');
[fail_h, fail_p] = ttest(fail_base, fail_resp, 'dim', 1, 'tail', 'both');
[press_h, press_p] = ttest(press_base, press_resp, 'dim', 1, 'tail', 'both');

%% 3. define cells by response to event
release_resp_cells = find(release_h);
success_resp_cells = find(success_h);
fail_resp_cells = find(fail_h);
press_resp_cells = find(press_h);
fail_only_cells = fail_resp_cells(find(ismember(fail_resp_cells, success_resp_cells)==0));
success_only_cells = success_resp_cells(find(ismember(success_resp_cells, fail_resp_cells)==0));
fail_and_success_cells = find(and(fail_h, success_h));
press_and_success_cells = find(and(success_h,press_h));
press_notsuccess_cells = press_resp_cells(find(ismember(press_resp_cells, success_resp_cells)==0));
success_notpress_cells = success_resp_cells(find(ismember(success_resp_cells, press_resp_cells)==0));
noresponse_cells = find((success_h+fail_h+press_h)==0);
save([dest_sub '_cell_categories.mat'], 'release_resp_cells', 'success_resp_cells', 'fail_resp_cells', 'press_resp_cells','noresponse_cells', 'fail_only_cells', 'success_only_cells', 'fail_and_success_cells', 'press_and_success_cells', 'press_notsuccess_cells', 'success_notpress_cells');

%test success vs fail responses (in cells that respond to both)
[fail_v_success_h, fail_v_success_p] = ttest(fail_resp_avg(1,fail_and_success_cells),success_resp_avg(1,fail_and_success_cells), 'tail', 'left');
save([dest_sub '_pvals.mat'], 'fail_v_success_p', 'release_p','success_p', 'fail_p', 'press_p','release_h','success_h', 'fail_h', 'press_h');

%% 4. plotting
tt =((-pre_release_frames:post_release_frames).*double(ifi))./1000;
%overlay average success/failure for each ROI
figure;
avg_all = [avg_success avg_fail];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(n,n2,ic)
    errorbar(tt,avg_success(ic,:), sem_success(ic,:),'k')
    hold on;
    errorbar(tt,avg_fail(ic,:), sem_fail(ic,:),'r')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle([date ' ' mouse ' Average release: Success- black (n = ' num2str(size(success_movie,1)) ' trials); Failure- red (n = ' num2str(size(fail_movie,1)) ' trials)'])
orient landscape
print([dest_sub '_success_fail_allTCs.eps'], '-depsc');
print([dest_sub '_success_fail_allTCs.pdf'], '-dpdf');
%overlay average success/failure across ROIs
figure; 
subplot(2,2,1)
errorbar(tt,avg_success_all, sem_success_all,'k')
hold on;
errorbar(tt,avg_fail_all, sem_fail_all,'r')
title(['Avg all cells; n = ' num2str(size(avg_success,1))])
ylabel('dF/F')
xlabel('Time (ms)')
subplot(2,2,2)
errorbar(tt,mean(avg_success(find(release_h),:),1),std(avg_success(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
hold on;
errorbar(tt,mean(avg_fail(find(release_h),:),1),std(avg_fail(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'r')
title(['Avg release resp p<0.05; n = ' num2str(sum(release_h,2))])
ylabel('dF/F')
xlabel('Time (ms)')
subplot(2,2,3)
errorbar(tt,mean(avg_success(find(success_h),:),1),std(avg_success(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'k')
hold on;
errorbar(tt,mean(avg_fail(find(success_h),:),1),std(avg_fail(find(success_h),:),[],1)./sqrt(sum(success_h,2)),'r')
title(['Avg success resp p<0.05; n = ' num2str(sum(success_h,2))])
ylabel('dF/F')
xlabel('Time (ms)')
suptitle([date ' ' mouse ' Average release: Success- black; Failure- red'])
print([dest_sub '_success_fail_avgTCs.eps'], '-depsc');
print([dest_sub '_success_fail_avgTCs.pdf'], '-dpdf');

%overlay average press/release for each ROI
figure;
avg_all = [avg_press avg_release];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(n,n2,ic)
    errorbar(tt,avg_press(ic,:), sem_press(ic,:),'c')
    hold on;
    errorbar(tt,avg_release(ic,:), sem_release(ic,:),'k')
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle([date ' ' mouse ' Average release: black (n = ' num2str(size(release_movie,1)) ' trials); Press- cyan (n = ' num2str(size(press_long_movie,1)) ' trials)'])
orient landscape
print([dest_sub '_press_release_allTCs.eps'], '-depsc');
print([dest_sub '_press_release_allTCs.pdf'], '-dpdf');
%overlay average press/release across ROIs
figure;
subplot(2,2,1)
errorbar(tt,avg_release_all, sem_release_all,'k')
hold on;
errorbar(tt,avg_press_all, sem_press_all,'c')
title(['Avg all cells; n = ' num2str(size(avg_release,1))])
ylabel('dF/F')
xlabel('Time (ms)')
subplot(2,2,2)
errorbar(tt,mean(avg_release(find(release_h),:),1),std(avg_release(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'k')
hold on;
errorbar(tt,mean(avg_press(find(release_h),:),1),std(avg_press(find(release_h),:),[],1)./sqrt(sum(release_h,2)),'c')
title(['Avg release resp p<0.05; n = ' num2str(sum(release_h,2))])
ylabel('dF/F')
xlabel('Time (ms)')
subplot(2,2,3)
errorbar(tt,mean(avg_release(find(press_h),:),1),std(avg_release(find(press_h),:),[],1)./sqrt(sum(press_h,2)),'k')
hold on;
errorbar(tt,mean(avg_press(find(press_h),:),1),std(avg_press(find(press_h),:),[],1)./sqrt(sum(press_h,2)),'c')
title(['Avg press resp p<0.05; n = ' num2str(sum(press_h,2))])
ylabel('dF/F')
xlabel('Time (ms)')
suptitle([date ' ' mouse ' Average release: Release- black; Press- cyan'])
print([dest_sub '_press_release_avgTCs.eps'], '-depsc');
print([dest_sub '_press_release_avgTCs.pdf'], '-dpdf');

figure;
all_resp = [success_resp_avg fail_resp_avg press_resp_avg];
cmax = max(all_resp,[],2);
cmin = min(all_resp,[],2);
subplot(2,2,1)
mask_final_temp = zeros(size(mask_final));
mask_final_temp(find(mask_final>0)) = -1;
imagesq(reshape(mask_final_temp, [sz(1) sz(2)]));
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
colorbar

success_mask = mask_final;
for ic = 1:nCells
    if success_h(ic)
        success_mask(find(success_mask==ic))=success_resp_avg(1,ic);
    else
        success_mask(find(success_mask==ic))=0;
    end
end
success_mask = reshape(success_mask,[sz(1) sz(2)]);
subplot(2,2,2)
imagesq(success_mask);
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
title('success')
colorbar

fail_mask = mask_final;
for ic = 1:nCells
    if fail_h(ic)
        fail_mask(find(fail_mask==ic))=fail_resp_avg(1,ic);
    else
        fail_mask(find(fail_mask==ic))=0;
    end
end
fail_mask = reshape(fail_mask,[sz(1) sz(2)]);
subplot(2,2,3)
imagesq(fail_mask);
clim([cmin cmax])
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
title('failure')
colorbar

press_mask = mask_final;
for ic = 1:nCells
    if press_h(ic)
        press_mask(find(press_mask==ic))=press_resp_avg(1,ic);
    else
        press_mask(find(press_mask==ic))=0;
    end
end
press_mask = reshape(press_mask,[sz(1) sz(2)]);
subplot(2,2,4)
imagesq(press_mask);
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
colorbar
title('press')

suptitle([mouse ' ' date ' cell responses']);
print([dest_sub '_cell_responses_FOV.eps'], '-depsc');
print([dest_sub '_cell_responses_FOV.pdf'], '-dpdf');

%plot amplitude response by cell
lever = strvcat('press', 'release','success', 'fail');
figure;
subplot(2,2,1)
for ic = 1:nCells
    plot(1:4, [press_resp_avg(1,ic) release_resp_avg(1,ic) success_resp_avg(1,ic) fail_resp_avg(1,ic)], '-ok')
    hold on
    if press_h(1,ic)
        plot(1,press_resp_avg(1,ic),'or')
        hold on
    end
    if release_h(1,ic)
        plot(2,release_resp_avg(1,ic),'or')
        hold on
    end
    if success_h(1,ic)
        plot(3,success_resp_avg(1,ic),'or')
        hold on
    end
    if fail_h(1,ic)
        plot(4,fail_resp_avg(1,ic),'or')
        hold on
    end
end
set(gca, 'XTick', 1:4, 'XTickLabel', lever);
xlim([0.5 4.5])
title('Resp amp- red: p<0.05')
ylabel('dF/F')

subplot(2,2,2)
errorbar(1:4, [mean(press_resp_avg,2) mean(release_resp_avg,2) mean(success_resp_avg,2) mean(fail_resp_avg,2)], [std(press_resp_avg,[],2) std(release_resp_avg,[],2) std(success_resp_avg,[],2) std(fail_resp_avg,[],2)]./sqrt(nCells),'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', lever);
xlim([0.5 4.5])
ylabel('dF/F')
title(['Avg resp- all cells: n = ' num2str(nCells)])
subplot(2,2,3)
errorbar(1:4, [mean(press_resp_avg(1,find(release_h)),2) mean(release_resp_avg(1,find(release_h)),2) mean(success_resp_avg(1,find(release_h)),2) mean(fail_resp_avg(1,find(release_h)),2)], [std(press_resp_avg(1,find(release_h)),[],2) std(release_resp_avg(1,find(release_h)),[],2) std(success_resp_avg(1,find(release_h)),[],2) std(fail_resp_avg(1,find(release_h)),[],2)]./sqrt(sum(release_h,2)),'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', lever);
xlim([0.5 4.5])
ylabel('dF/F')
title(['Avg resp- release resp cells: n = ' num2str(sum(release_h,2))])
subplot(2,2,4)
errorbar(1:4, [mean(press_resp_avg(1,find(press_h)),2) mean(release_resp_avg(1,find(press_h)),2) mean(success_resp_avg(1,find(press_h)),2) mean(fail_resp_avg(1,find(press_h)),2)], [std(press_resp_avg(1,find(press_h)),[],2) std(release_resp_avg(1,find(press_h)),[],2) std(success_resp_avg(1,find(press_h)),[],2) std(fail_resp_avg(1,find(press_h)),[],2)]./sqrt(sum(press_h,2)),'ok')
set(gca, 'XTick', 1:4, 'XTickLabel', lever);
xlim([0.5 4.5])
ylabel('dF/F')
title(['Avg resp- press resp cells: n = ' num2str(sum(press_h,2))])
suptitle([date ' ' mouse 'avg resp amplitude'])
print([dest_sub '_cell_resp_amplitude.eps'], '-depsc');
print([dest_sub '_cell_resp_amplitude.pdf'], '-dpdf');
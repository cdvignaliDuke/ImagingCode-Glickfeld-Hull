% finds and maps responsive cells for each condition
% 1. calculate average timecourses for NR/release events
% 2. calculate variability by trial over base and resp windows
% 3. ttest for significant responses
% 4. define cells by response to event

% load([dest 'ROI_TCs.mat']);
% load([dest '_release_movies.mat'])
% load([dest '_cue_movies.mat'])
nCells = size(NR_movie,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n) %function
    n2= n-1;
else
    n2 = n;
end

%% 1. calculate average timecourses for NR/release events

%average and sem across NRes
avg_NR = squeeze(mean(NR_movie,1));
sem_NR = squeeze(std(NR_movie,1)./sqrt(size(NR_movie,1)));
%average and sem across ROIs (and trials)
avg_NR_all = mean(avg_NR,1);
sem_NR_all = std(avg_NR,1)./sqrt(size(avg_NR,1));

%average and sem across NRes
avg_OR = squeeze(mean(OR_movie,1));
sem_OR = squeeze(std(OR_movie,1)./sqrt(size(OR_movie,1)));
%average and sem across ROIs (and trials)
avg_OR_all = mean(avg_OR,1);
sem_OR_all = std(avg_OR,1)./sqrt(size(avg_OR,1));

%average and sem across NRes
avg_UR = squeeze(mean(UR_movie,1));
sem_UR = squeeze(std(UR_movie,1)./sqrt(size(UR_movie,1)));
%average and sem across ROIs (and trials)
avg_UR_all = mean(avg_UR,1);
sem_UR_all = std(avg_UR,1)./sqrt(size(avg_UR,1));

save([dest '_cell_TCs.mat'], 'avg_NR', 'sem_NR', 'avg_OR', 'sem_OR', 'avg_UR', 'sem_UR')
%% 2. calculate response amplitude and variablity by trial over base and resp windows

%by NR
[NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem] = findRespCell(NR_movie, pre_cue_frames, ifi, dest);

[OR_h, OR_p, OR_resp_cells, OR_resp_avg, OR_resp_sem] = findRespCell(OR_movie, pre_cue_frames, ifi, dest);

[UR_h, UR_p, UR_resp_cells, UR_resp_avg, UR_resp_sem] = findRespCell(UR_movie, pre_cue_frames, ifi, dest);

noresponse_cells = find((NR_h+OR_h+UR_h)==0);

if ~isempty(OR_movie)
    
    allresp_cells = NR_h & OR_h;

    NOR_h = allresp_cells;
else
    allresp_cells = NR_h;
end

NR_h = logical(NR_h);

if ~isempty(OR_movie)
   OR_h = logical(OR_h);
else
   OR_h = zeros(size(OR_h));
end

if ~isempty(UR_movie)
    UR_h = logical(UR_h);
    allresp_cells = allresp_cells & UR_h;
    NUR_h = NR_h & UR_h;
    OUR_h = OR_h & UR_h;
end

save([dest '_cell_categories.mat'], 'NR_resp_cells', 'OR_resp_cells', 'UR_resp_cells', 'noresponse_cells', 'allresp_cells');

save([dest '_pvals.mat'], 'NR_p', 'NR_h', 'OR_p', 'OR_h', 'UR_p', 'UR_h');


%% 4. plotting
tt =((-pre_cue_frames:post_cue_frames).*double(ifi))./1000;
%overlay average success/failure for each ROI
figure;
avg_all = [avg_NR avg_OR avg_UR];
ymax = max(max(avg_all,[],2),[],1);
ymin = min(min(avg_all,[],2),[],1);
for ic = 1:nCells
    subplot(n,n2,ic)
%     errorbar(tt,avg_NR(ic,:), sem_NR(ic,:),'k');hold on;
%     errorbar(tt,avg_OR(ic,:), sem_OR(ic,:),'r');
%     errorbar(tt,avg_UR(ic,:), sem_UR(ic,:),'b');
    plot(tt,avg_NR(ic,:),'k');hold on;
    if ~isempty(OR_movie)
        plot(tt,avg_OR(ic,:),'r');
    end
    if ~isempty(UR_movie)
        plot(tt,avg_UR(ic,:),'g');
    end
    ylim([ymin*1.1 ymax*1.1])
    xlim([tt(1) tt(end)])
end
suptitle([date ' ' mouse ' Average cue resp: Normal- black (n = ' num2str(size(NR_movie,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie,1)) ' trials)'])
orient landscape
print([dest '_allTCs.eps'], '-depsc');
print([dest '_allTCs.pdf'], '-dpdf');

%overlay average success/failure across ROIs
h=figure; 
subplot(3,3,1)
errorbar(tt,avg_NR_all, sem_NR_all,'k')
hold on;
if ~isempty(OR_movie)
    errorbar(tt,avg_OR_all, sem_OR_all,'r');
end
if ~isempty(UR_movie)
    errorbar(tt,avg_UR_all, sem_UR_all,'g')
end
% bar(tt,mean(lick_trace_NR)/5, 'k');
% bar(tt,mean(lick_trace_OR)/5, 'r');
% if ~isempty(UR_movie)
%     bar(tt,mean(lick_trace_UR), 'g-');
% end
title(['Avg all cells; n = ' num2str(size(avg_NR,1))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

% plot NR resp cells
subplot(3,3,2)
errorbar(tt,mean(avg_NR(NR_h,:),1),std(avg_NR(NR_h,:),[],1)./sqrt(sum(NR_h,2)),'k')
hold on;

bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');

title(['Avg Normal Reward resp cells p<0.05; n = ' num2str(sum(NR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

% plot OR resp cells
if ~isempty(OR_movie)
subplot(3,3,3)
errorbar(tt,mean(avg_OR(OR_h,:),1),std(avg_OR(OR_h,:),[],1)./sqrt(sum(OR_h,2)),'r')
hold on;

bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');

title(['Avg Omitted Reward resp cells p<0.05; n = ' num2str(sum(OR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot UR resp cells
if ~isempty(UR_movie)
subplot(3,3,4)
errorbar(tt,mean(avg_UR(UR_h,:),1),std(avg_UR(UR_h,:),[],1)./sqrt(sum(UR_h,2)),'g')
hold on;

bar(tt,mean(lick_trace_UR)/8, 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');

title(['Avg Unexpected Reward resp cells p<0.05; n = ' num2str(sum(UR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot NR and OR resp cells
if ~isempty(OR_movie)
subplot(3,3,5)
errorbar(tt,mean(avg_NR(NOR_h,:),1),std(avg_NR(NOR_h,:),[],1)./sqrt(sum(NOR_h,2)),'k')
hold on;
errorbar(tt,mean(avg_OR(NOR_h,:),1),std(avg_OR(NOR_h,:),[],1)./sqrt(sum(NOR_h,2)),'r')
% bar(tt,mean(lick_trace_NR)/8, 'k');

title(['Avg Normal and Omitted Reward resp cells p<0.05; n = ' num2str(sum(NOR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot NR and UR resp cells
if ~isempty(UR_movie)
subplot(3,3,6)
errorbar(tt,mean(avg_NR(NUR_h,:),1),std(avg_NR(NUR_h,:),[],1)./sqrt(sum(NUR_h,2)),'k')
hold on;
errorbar(tt,mean(avg_UR(NUR_h,:),1),std(avg_UR(NUR_h,:),[],1)./sqrt(sum(NUR_h,2)),'g')
% bar(tt,mean(lick_trace_NR)/8, 'k');

title(['Avg Normal and Unexpected Reward resp cells p<0.05; n = ' num2str(sum(NUR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

% plot OR and UR resp cells
if ~isempty(UR_movie) && ~isempty(OR_movie)
subplot(3,3,7)
errorbar(tt,mean(avg_OR(OUR_h,:),1),std(avg_OR(OUR_h,:),[],1)./sqrt(sum(OUR_h,2)),'r')
hold on;
errorbar(tt,mean(avg_UR(OUR_h,:),1),std(avg_UR(OUR_h,:),[],1)./sqrt(sum(OUR_h,2)),'g')
% bar(tt,mean(lick_trace_NR)/8, 'k');

title(['Avg Omitted and Unexpected Reward resp cells p<0.05; n = ' num2str(sum(OUR_h,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);
end

subplot(3,3,8)
errorbar(tt,mean(avg_NR(find(allresp_cells),:),1),std(avg_NR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'k')
hold on;
errorbar(tt,mean(avg_OR(find(allresp_cells),:),1),std(avg_OR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'r')
errorbar(tt,mean(avg_UR(find(allresp_cells),:),1),std(avg_UR(find(allresp_cells),:),[],1)./sqrt(sum(allresp_cells,2)),'g')

bar(tt,mean(lick_trace_NR)/8, 'k', 'ShowBaseLine', 'off');
if ~isempty(OR_movie)
bar(tt,mean(lick_trace_OR)/8, 'FaceColor', [1,0,0], 'EdgeColor', [1,0,0], 'ShowBaseLine', 'off');
end

if ~isempty(UR_movie)
    bar(tt,mean(lick_trace_UR/8), 'FaceColor', [0,1,0], 'EdgeColor', [0,1,0], 'ShowBaseLine', 'off');
end

title(['Avg all resp cells p<0.05; n = ' num2str(sum(allresp_cells,2))])
ylabel('dF/F')
xlabel('Time (s)')
xlim([-floor(pre_cue_frames/ifi) floor(post_cue_frames/ifi)]);

suptitle([date ' ' mouse ' Average cue resp: Normal- black; Omit- red; Unexpect- green bars- licking data'])
print([dest '_avgTCs.eps'], '-depsc');
print([dest '_avgTCs.pdf'], '-dpdf');
savefig(h, [dest '_avgTCs']); 

figure;
if ~isempty(UR_movie) && ~isempty(OR_movie)
    all_resp = [NR_resp_avg OR_resp_avg UR_resp_avg];
elseif isempty(UR_movie) && ~isempty(OR_movie)
    all_resp = [NR_resp_avg OR_resp_avg];
elseif ~isempty(UR_movie) && isempty(OR_movie)
    all_resp = [NR_resp_avg UR_resp_avg];
elseif isempty(UR_movie) && isempty(OR_movie)
    all_resp = [NR_resp_avg];
end
cmax = max(all_resp,[],2);
cmin = min(all_resp,[],2);
subplot(2,2,1)
mask_final_temp = zeros(size(mask_final));
mask_final_temp(find(mask_final>0)) = 1;
imagesq(reshape(mask_final_temp, [sz(1) sz(2)]));
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
colorbar

NR_mask = mask_final;
for ic = 1:nCells
    if NR_h(ic)
        NR_mask(find(NR_mask==ic))=NR_resp_avg(1,ic);
    else
        NR_mask(find(NR_mask==ic))=0;
    end
end
NR_mask = reshape(NR_mask,[sz(1) sz(2)]);
subplot(2,2,2)
imagesq(NR_mask);
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
title('normal')
colorbar

if ~isempty(OR_movie)
OR_mask = mask_final;
for ic = 1:nCells
    if OR_h(ic)
        OR_mask(find(OR_mask==ic))=OR_resp_avg(1,ic);
    else
        OR_mask(find(OR_mask==ic))=0;
    end
end
OR_mask = reshape(OR_mask,[sz(1) sz(2)]);
subplot(2,2,3)
imagesq(OR_mask);
clim([cmin cmax])
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
title('Omit')
colorbar
end

if ~isempty(UR_movie)
UR_mask = mask_final;
for ic = 1:nCells
    if UR_h(ic)
        UR_mask(find(UR_mask==ic))=UR_resp_avg(1,ic);
    else
        UR_mask(find(UR_mask==ic))=0;
    end
end
UR_mask = reshape(UR_mask,[sz(1) sz(2)]);
subplot(2,2,4)
imagesq(UR_mask);
clim([cmin cmax]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
colorbar
title('Unexpect')
end

suptitle([mouse ' ' date ' cell responses']);
print([dest '_cell_responses_FOV.eps'], '-depsc');
print([dest '_cell_responses_FOV.pdf'], '-dpdf');

%plot amplitude response by cell
% lever = strvcat('NR', 'release','success', 'fail');
% figure;
% subplot(2,2,1)
% for ic = 1:nCells
%     plot(1:4, [NR_resp_avg(1,ic) release_resp_avg(1,ic) success_resp_avg(1,ic) fail_resp_avg(1,ic)], '-ok')
%     hold on
%     if NR_h(1,ic)
%         plot(1,NR_resp_avg(1,ic),'or')
%         hold on
%     end
%     if release_h(1,ic)
%         plot(2,release_resp_avg(1,ic),'or')
%         hold on
%     end
%     if success_h(1,ic)
%         plot(3,success_resp_avg(1,ic),'or')
%         hold on
%     end
%     if fail_h(1,ic)
%         plot(4,fail_resp_avg(1,ic),'or')
%         hold on
%     end
% end
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% title('Resp amp- red: p<0.05')
% ylabel('dF/F')
% 
% subplot(2,2,2)
% errorbar(1:4, [mean(NR_resp_avg,2) mean(release_resp_avg,2) mean(success_resp_avg,2) mean(fail_resp_avg,2)], [std(NR_resp_avg,[],2) std(release_resp_avg,[],2) std(success_resp_avg,[],2) std(fail_resp_avg,[],2)]./sqrt(nCells),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- all cells: n = ' num2str(nCells)])
% subplot(2,2,3)
% errorbar(1:4, [mean(NR_resp_avg(1,find(release_h)),2) mean(release_resp_avg(1,find(release_h)),2) mean(success_resp_avg(1,find(release_h)),2) mean(fail_resp_avg(1,find(release_h)),2)], [std(NR_resp_avg(1,find(release_h)),[],2) std(release_resp_avg(1,find(release_h)),[],2) std(success_resp_avg(1,find(release_h)),[],2) std(fail_resp_avg(1,find(release_h)),[],2)]./sqrt(sum(release_h,2)),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- release resp cells: n = ' num2str(sum(release_h,2))])
% subplot(2,2,4)
% errorbar(1:4, [mean(NR_resp_avg(1,find(NR_h)),2) mean(release_resp_avg(1,find(NR_h)),2) mean(success_resp_avg(1,find(NR_h)),2) mean(fail_resp_avg(1,find(NR_h)),2)], [std(NR_resp_avg(1,find(NR_h)),[],2) std(release_resp_avg(1,find(NR_h)),[],2) std(success_resp_avg(1,find(NR_h)),[],2) std(fail_resp_avg(1,find(NR_h)),[],2)]./sqrt(sum(NR_h,2)),'ok')
% set(gca, 'XTick', 1:4, 'XTickLabel', lever);
% xlim([0.5 4.5])
% ylabel('dF/F')
% title(['Avg resp- NR resp cells: n = ' num2str(sum(NR_h,2))])
% suptitle([date ' ' mouse 'avg resp amplitude'])
% print([dest '_cell_resp_amplitude.eps'], '-depsc');
% print([dest '_cell_resp_amplitude.pdf'], '-dpdf');



clear all
close all
CRP_expt_list
nexp = size(expt(2).date,1);
rew_all_d1 = [];
rew_all_d2 = [];
omit_all_d1 = [];
omit_all_d2 = [];
mouse_use = [];
dend_num = [];
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
share_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\ClassicConditioningPaper';
for iexp = 1:nexp
    id = 2;
    mouse = expt(id).mouse(iexp,:);
    date = expt(id).date(iexp,:);
    run = expt(id).run(iexp,:);
    fprintf([mouse ' ' date '\n'])
    img_fn = [date '_' mouse];
    if exist(fullfile(lg_out, img_fn, [img_fn '_D2toD1_overlap.mat']))
        load(fullfile(lg_out, img_fn, [img_fn '_D2toD1_overlap.mat']))
        if sum(~isnan(overlap_id))
            load(fullfile(lg_out, img_fn, [img_fn '_targetAlign.mat']))
            ind = find(~isnan(overlap_id));
            rew_all_d2 = [rew_all_d2 nanmean(targetAlign_events(:,ind,ind_rew),3)];
            omit_all_d2 = [omit_all_d2 nanmean(targetAlign_events(:,ind,ind_omit),3)];
            id = 1;
            nexp1 = size(expt(1).date,1);
            x = 0;
            for iexp1 = 1:nexp1
                if strcmp(expt(id).mouse(iexp1,:),mouse)
                    x = 1;
                    break
                end
            end
            if x == 1
                date = expt(id).date(iexp1,:);
                run = expt(id).run(iexp1,:);
                fprintf([mouse ' ' date '\n'])
                img_fn = [date '_' mouse];
                load(fullfile(lg_out, img_fn, [img_fn '_targetAlign.mat']))
                ind1 = overlap_id(ind);
                rew_all_d1 = [rew_all_d1 nanmean(targetAlign_events(:,ind1,ind_rew),3)];
                omit_all_d1 = [omit_all_d1 nanmean(targetAlign_events(:,ind1,ind_omit),3)];
                mouse_use = [mouse_use; mouse];
                dend_num = [dend_num; length(ind1)];
            end
        end
    end
end
save([share_out '\CC_crossday\D1vD2_matchedCell_Summary.mat'], 'rew_all_d1','omit_all_d1', 'rew_all_d2','omit_all_d2', 'mouse_use')

rewdelay_frames = round(0.6.*frameRateHz);
pre_rew_win = cell(1,2);
post_rew_win = cell(1,2);
for id = 1:2
    load(fullfile(lg_out,'CC_Summary',['Day' num2str(id)],['Day' num2str(id) '_PeakSpikingResp_Windows.mat']))
    pre_rew_win{id} = preresp_rew_range;
    post_rew_win{id} = postresp_rew_range;
end
rew_all_d1_sub = rew_all_d1-mean(rew_all_d1(1:prewin_frames,:),1);
rew_all_d2_sub = rew_all_d2-mean(rew_all_d2(1:prewin_frames,:),1);

avg_d1_rew = mean(rew_all_d1,2);
avg_d2_rew = mean(rew_all_d2,2);
base_d1_avg = mean(mean(rew_all_d1(ceil(prewin_frames./2):prewin_frames,:),1),2);
base_d2_avg = mean(mean(rew_all_d2(ceil(prewin_frames./2):prewin_frames,:),1),2);
base_d1_std = std(mean(rew_all_d1(ceil(prewin_frames./2):prewin_frames,:),2),[],1);
base_d2_std = std(mean(rew_all_d2(ceil(prewin_frames./2):prewin_frames,:),2),[],1);
[prerew_max_rew prerew_max_ind_rew] = max(avg_d2_rew(prewin_frames+1:prewin_frames+rewdelay_frames+1,:),[],1);
[postrew_max_rew postrew_max_ind_rew] = max(avg_d1_rew(prewin_frames+rewdelay_frames+1:prewin_frames+rewdelay_frames+1+rewdelay_frames,:),[],1);
[base_max_d1 base_max_ind_d1] = max(avg_d1_rew(prewin_frames-rewdelay_frames-rewdelay_frames:prewin_frames,:),[],1);
[base_max_d2 base_max_ind_d2] = max(avg_d2_rew(prewin_frames-rewdelay_frames-rewdelay_frames:prewin_frames,:),[],1);
 
nIC = size(rew_all_d1,2);
figure;
subplot(4,1,1)
shadedErrorBar(tt, nanmean(rew_all_d1,2).*1000/frameRateHz, (nanstd(rew_all_d1,[],2)./sqrt(nIC)).*1000/frameRateHz,'k');
hold on
shadedErrorBar(tt, nanmean(rew_all_d2,2).*1000/frameRateHz, (nanstd(rew_all_d2,[],2)./sqrt(nIC)).*1000/frameRateHz,'b');
ylabel('Spike rate (Hz)')
xlabel('Time from cue (ms)')
xlim([-500 2000])
ylim([0 3])
vline([tt(pre_rew_win{2}{3}(2)) tt(post_rew_win{1}{3}(2))],'k')
if prerew_max_rew > base_d2_avg+(base_d2_std.*4)
	vline(tt(prewin_frames+prerew_max_ind_rew),'r')
else
    vline(tt(prewin_frames+prerew_max_ind_rew),'--r')
end
if postrew_max_rew > base_d1_avg+(base_d1_std.*4)
	vline(tt(prewin_frames+rewdelay_frames+postrew_max_ind_rew),'r')
else
    vline(tt(prewin_frames+rewdelay_frames+postrew_max_ind_rew),'--r')
end
title('Reward')
subplot(4,2,3)
scatter(mean(rew_all_d1_sub(pre_rew_win{2}{3},:),1).*1000/frameRateHz,mean(rew_all_d2_sub(pre_rew_win{2}{3},:),1).*1000/frameRateHz,'ok')
hold on
errorbarxy(mean(mean(rew_all_d1_sub(pre_rew_win{2}{3},:),1),2).*1000/frameRateHz,mean(mean(rew_all_d2_sub(pre_rew_win{2}{3},:),1),2).*1000/frameRateHz,std(mean(rew_all_d1_sub(pre_rew_win{2}{3},:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),std(mean(rew_all_d2_sub(pre_rew_win{2}{3},:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),{'or','r','r'})
xlabel('Peak Firing Rate - D1')
ylabel('Peak Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre] = ttest(mean(rew_all_d1_sub(pre_rew_win{2}{3},:),1).*1000/frameRateHz,mean(rew_all_d2_sub(pre_rew_win{2}{3},:),1).*1000/frameRateHz,'tail','both');
title(['Pre reward - D2 peak- p = ' num2str(chop(p_pre,2))])
subplot(4,2,4)
scatter(mean(rew_all_d1_sub(post_rew_win{1}{3},:),1).*1000/frameRateHz,mean(rew_all_d2_sub(post_rew_win{1}{3},:),1).*1000/frameRateHz,'ok')
hold on
errorbarxy(mean(mean(rew_all_d1_sub(post_rew_win{1}{3},:),1),2).*1000/frameRateHz,mean(mean(rew_all_d2_sub(post_rew_win{1}{3},:),1),2).*1000/frameRateHz,std(mean(rew_all_d1_sub(post_rew_win{1}{3},:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),std(mean(rew_all_d2_sub(post_rew_win{1}{3},:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),{'or','r','r'})
xlabel('Peak Firing Rate - D1')
ylabel('Peak Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_post p_post] = ttest(mean(rew_all_d1_sub(post_rew_win{1}{3},:),1).*1000/frameRateHz,mean(rew_all_d2_sub(post_rew_win{1}{3},:),1).*1000/frameRateHz);
title(['Post reward - D1 peak- p = ' num2str(chop(p_post,2))])

preresp_rew_range = [prewin_frames+prerew_max_ind_rew-1:prewin_frames+prerew_max_ind_rew+1]; 
postresp_rew_range = [prewin_frames+rewdelay_frames+postrew_max_ind_rew-1:prewin_frames+rewdelay_frames+postrew_max_ind_rew+1];
subplot(4,2,5)
scatter(mean(rew_all_d1_sub(preresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(preresp_rew_range,:),1).*1000/frameRateHz,'ok')
hold on
errorbarxy(mean(mean(rew_all_d1_sub(preresp_rew_range,:),1),2).*1000/frameRateHz,mean(mean(rew_all_d2_sub(preresp_rew_range,:),1),2).*1000/frameRateHz,std(mean(rew_all_d1_sub(preresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),std(mean(rew_all_d2_sub(preresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),{'or','r','r'})
xlabel('Peak Firing Rate - D1')
ylabel('Peak Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre2] = ttest(mean(rew_all_d1_sub(preresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(preresp_rew_range,:),1).*1000/frameRateHz);
title(['Pre reward - D1 peak- p = ' num2str(chop(p_pre2,2))])
subplot(4,2,6)
scatter(mean(rew_all_d1_sub(postresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(postresp_rew_range,:),1).*1000/frameRateHz,'ok')
hold on
errorbarxy(mean(mean(rew_all_d1_sub(postresp_rew_range,:),1),2).*1000/frameRateHz,mean(mean(rew_all_d2_sub(postresp_rew_range,:),1),2).*1000/frameRateHz,std(mean(rew_all_d1_sub(postresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),std(mean(rew_all_d2_sub(postresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),{'or','r','r'})
xlabel('Peak Firing Rate - D1')
ylabel('Peak Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_post p_post2] = ttest(mean(rew_all_d1_sub(postresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(postresp_rew_range,:),1).*1000/frameRateHz);
title(['Post reward - D2 peak- p = ' num2str(chop(p_post2,2))])
subplot(4,2,7)
scatter(mean(rew_all_d1_sub(postresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(preresp_rew_range,:),1).*1000/frameRateHz,'ok')
hold on
errorbarxy(mean(mean(rew_all_d1_sub(postresp_rew_range,:),1),2).*1000/frameRateHz,mean(mean(rew_all_d2_sub(preresp_rew_range,:),1),2).*1000/frameRateHz,std(mean(rew_all_d1_sub(postresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),std(mean(rew_all_d2_sub(postresp_rew_range,:),1),[],2).*1000/frameRateHz./sqrt(size(rew_all_d1,2)),{'or','r','r'})
xlabel('D1- Post reward')
ylabel('D2- Pre reward')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_post p_post2] = ttest(mean(rew_all_d1_sub(postresp_rew_range,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(preresp_rew_range,:),1).*1000/frameRateHz);
title(['p = ' num2str(chop(p_post2,2))])

suptitle(['Matched cells (n = ' num2str(nIC) ') - D1 (black) vs D2 (blue)'])
savefig([share_out '\CC_crossday\D1vD2_matchedCell_Summary_new.fig'])

figure;
subplot(2,2,1)
scatter(mean(rew_all_d1_sub(prewin_frames+prerew_max_ind_rew-1:prewin_frames+prerew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d1_sub(base_max_ind_d1-1:base_max_ind_d1+1,:),1).*1000/frameRateHz,'ok');
xlabel('Peak Firing Rate - D1')
ylabel('Base Firing Rate - D1')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre] = ttest(mean(rew_all_d1_sub(prewin_frames+prerew_max_ind_rew-1:prewin_frames+prerew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'tail','both');
title(['Pre reward - D1 peak- p = ' num2str(chop(p_pre,2))])
subplot(2,2,2)
scatter(mean(rew_all_d2_sub(prewin_frames+prerew_max_ind_rew-1:prewin_frames+prerew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'ok');
xlabel('Peak Firing Rate - D2')
ylabel('Base Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre] = ttest(mean(rew_all_d2_sub(prewin_frames+prerew_max_ind_rew-1:prewin_frames+prerew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'tail','both');
title(['Pre reward - D2 peak- p = ' num2str(chop(p_pre,2))])
subplot(2,2,3)
scatter(mean(rew_all_d1_sub(prewin_frames+rewdelay_frames+postrew_max_ind_rew-1:prewin_frames+rewdelay_frames+postrew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d1_sub(base_max_ind_d1-1:base_max_ind_d1+1,:),1).*1000/frameRateHz,'ok');
xlabel('Peak Firing Rate - D1')
ylabel('Base Firing Rate - D1')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre] = ttest(mean(rew_all_d1_sub(prewin_frames+rewdelay_frames+postrew_max_ind_rew-1:prewin_frames+rewdelay_frames+postrew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'tail','both');
title(['Post reward - D1 peak- p = ' num2str(chop(p_pre,2))])
subplot(2,2,4)
scatter(mean(rew_all_d2_sub(prewin_frames+rewdelay_frames+postrew_max_ind_rew-1:prewin_frames+rewdelay_frames+postrew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'ok');
xlabel('Peak Firing Rate - D2')
ylabel('Base Firing Rate - D2')
axis square
ylim([-1 4])
xlim([-1 4])
refline(1,0)
[h_pre p_pre] = ttest(mean(rew_all_d2_sub(prewin_frames+rewdelay_frames+postrew_max_ind_rew-1:prewin_frames+rewdelay_frames+postrew_max_ind_rew+1,:),1).*1000/frameRateHz,mean(rew_all_d2_sub(base_max_ind_d2-1:base_max_ind_d2+1,:),1).*1000/frameRateHz,'tail','both');
title(['Post reward - D2 peak- p = ' num2str(chop(p_pre,2))])
suptitle(['Matched cells (n = ' num2str(nIC) ') - Resp vs Base'])
print([share_out '\CC_crossday\D1vD2_matchedCell_Summary_BaseComp.pdf'],'-bestfit','-dpdf')
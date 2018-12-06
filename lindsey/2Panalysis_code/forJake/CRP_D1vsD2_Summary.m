clear all
close all
CRP_expt_list
nexp = size(expt(2).date,1);
rew_all_d1 = [];
rew_all_d2 = [];
omit_all_d1 = [];
omit_all_d2 = [];
rew_all_df_d1 = [];
rew_all_df_d2 = [];
omit_all_df_d1 = [];
omit_all_df_d2 = [];
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
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
            rew_all_df_d2 = [rew_all_df_d2 nanmean(targetAligndFoverF(:,ind,ind_rew),3)];
            omit_all_df_d2 = [omit_all_df_d2 nanmean(targetAligndFoverF(:,ind,ind_omit),3)];
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
                rew_all_df_d1 = [rew_all_df_d1 nanmean(targetAligndFoverF(:,ind1,ind_rew),3)];
                omit_all_df_d1 = [omit_all_df_d1 nanmean(targetAligndFoverF(:,ind1,ind_omit),3)];
            end
        end
    end
end
nIC = size(rew_all_d1,2);
figure;
subplot(2,2,1)
shadedErrorBar(tt, nanmean(rew_all_df_d1,2), (nanstd(rew_all_df_d1,[],2)./sqrt(nIC)),'k');
hold on
shadedErrorBar(tt, nanmean(rew_all_df_d2,2), (nanstd(rew_all_df_d2,[],2)./sqrt(nIC)),'b');
ylabel('dF/F')
xlabel('Time from cue (ms)')
xlim([-500 2000])
ylim([-0.01 0.1])
vline([600], 'b')
title('Reward')
subplot(2,2,2)
shadedErrorBar(tt, nanmean(rew_all_d1,2).*1000/frameRateHz, (nanstd(rew_all_d1,[],2)./sqrt(nIC)).*1000/frameRateHz,'k');
hold on
shadedErrorBar(tt, nanmean(rew_all_d2,2).*1000/frameRateHz, (nanstd(rew_all_d2,[],2)./sqrt(nIC)).*1000/frameRateHz,'b');
ylabel('Spike rate (Hz)')
xlabel('Time from cue (ms)')
xlim([-500 2000])
ylim([0 3])
vline([600], 'b')
title('Reward')
subplot(2,2,3)
shadedErrorBar(tt, nanmean(omit_all_df_d1,2), (nanstd(omit_all_df_d1,[],2)./sqrt(nIC)),'k');
hold on
shadedErrorBar(tt, nanmean(omit_all_df_d2,2), (nanstd(omit_all_df_d2,[],2)./sqrt(nIC)),'b');
ylabel('dF/F')
xlabel('Time from cue (ms)')
xlim([-500 2000])
ylim([-0.01 0.1])
vline([600], 'b')
title('Omission')
subplot(2,2,4)
shadedErrorBar(tt, nanmean(omit_all_d1,2).*1000/frameRateHz, (nanstd(omit_all_d1,[],2)./sqrt(nIC)).*1000/frameRateHz,'k');
hold on
shadedErrorBar(tt, nanmean(omit_all_d2,2).*1000/frameRateHz, (nanstd(omit_all_d2,[],2)./sqrt(nIC)).*1000/frameRateHz,'b');
ylabel('Spike rate (Hz)')
xlabel('Time from cue (ms)')
xlim([-500 2000])
ylim([0 3])
vline([600], 'b')
title('Omission')
suptitle(['Matched cells (n = ' num2str(nIC) ') - D1 (black) vs D2 (blue)'])
print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' 'D1vD2_matchedCell_Summary.pdf'],'-bestfit','-dpdf')
clear all
close all
frameRateHz = 30;
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake\CC_summary';
share_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\ClassicConditioningPaper';
for id= 1:2
    load(fullfile(lg_out,['Day' num2str(id)], ['Day' num2str(id) '_allArea_Data.mat']))
    eval(['all_rew_d' num2str(id) ' = all_rew;'])
    eval(['all_omit_d' num2str(id) ' = all_omit;'])
    eval(['all_area_id_d' num2str(id) ' = all_area_id;'])
    eval(['preresp_rew_range_d' num2str(id) ' = preresp_rew_range;'])
    eval(['postresp_rew_range_d' num2str(id) ' = postresp_rew_range;'])
end

area_list = {'C1','C2','LS'};
figure;
for i = 1:length(area_list)
    ind_d1 = find(all_area_id_d1 == i);
    ind_d2 = find(all_area_id_d2 == i);

    preRew_d1_avg = mean(all_rew_d1(preresp_rew_range_d1{i},ind_d1),1)- mean(all_rew_d1(1:prewin_frames,ind_d1),1);
    preRew_d2_avg = mean(all_rew_d2(preresp_rew_range_d2{i},ind_d2),1)- mean(all_rew_d2(1:prewin_frames,ind_d2),1);
    [h_LS_preRew, p_LS_preRew] = ttest2(preRew_d1_avg,preRew_d2_avg);
    subplot(1,length(area_list),i)
    scatter(ones(size(preRew_d1_avg)), preRew_d1_avg .* 1000/frameRateHz)
    hold on
    scatter(2.*ones(size(preRew_d2_avg)), preRew_d2_avg .* 1000/frameRateHz)
    title(['Area ' area_list(i) '- p = ' num2str(chop(p_LS_preRew,2))])
    xlabel('Day')
    ylabel('FR (Hz)')
    ylim([-2 10])
    xlim([0 3])
    axis square
end
suptitle('D1 vs D2 - pre reward response')
savefig([share_out '\CC_summary\D1vsD2\D1vsD2_Area' area_list{i} '_Summary_preRewResp.fig'])

figure;
for i = 1:length(area_list)
    ind_d1 = find(all_area_id_d1 == i);
    ind_d2 = find(all_area_id_d2 == i);

    postRew_d1_avg = mean(all_rew_d1(postresp_rew_range_d1{i},ind_d1),1)- mean(all_rew_d1(1:prewin_frames,ind_d1),1);
    postRew_d2_avg = mean(all_rew_d2(postresp_rew_range_d2{i},ind_d2),1)- mean(all_rew_d2(1:prewin_frames,ind_d2),1);
    [h_LS_postRew, p_LS_postRew] = ttest2(postRew_d1_avg,postRew_d2_avg);
    subplot(1,length(area_list),i)
    scatter(ones(size(postRew_d1_avg)), postRew_d1_avg .* 1000/frameRateHz)
    hold on
    scatter(2.*ones(size(postRew_d2_avg)), postRew_d2_avg .* 1000/frameRateHz)
    title(['Area ' area_list(i) '- p = ' num2str(chop(p_LS_postRew,2))])
    xlabel('Day')
    ylabel('FR (Hz)')
    ylim([-2 10])
    xlim([0 3])
    axis square
end
suptitle('D1 vs D2 - post reward response')
savefig([share_out '\CC_summary\D1vsD2\D1vsD2_Area' area_list{i} '_Summary_postRewResp.fig'])


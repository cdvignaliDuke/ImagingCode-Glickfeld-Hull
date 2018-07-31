mouse_mat = strvcat('i698', 'i699','i671');
date_mat = strvcat('170923', '171013','170127');
run_str_mat = {'runs-001-003','runs-001','runs-002-004'};
frame_rate  =30;
nexp = size(mouse_mat,1);
mouse_str = [];
for iexp = 1:nexp
    mouse_str = [mouse_str mouse_mat(iexp,:)];
end

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% collect data
rad_targ_miss_all = zeros(1, size(mouse_mat,1));
rad_targ_hit_all = zeros(1, size(mouse_mat,1));
rad_targ_isi_all = zeros(3, size(mouse_mat,1));
for imouse = 1:size(mouse_mat,1)
    mouse = mouse_mat(imouse,:);
    date = date_mat(imouse,:);
    run_str = run_str_mat{imouse};
    
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eyeData.mat']))
    
    if imouse == 2
        run_str = 'runs-001-002';
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        nTrials = 225;
        MIx = MIx(1:nTrials);
        SIx = SIx(1:nTrials);        
        run_str = 'runs-001';
    else
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    end
    
    targetFramesOff = nan(1,nTrials);
    for itrial = 1:nTrials
        if MIx(itrial) || SIx(itrial)
            targetFramesOff(1,itrial) = tFramesOff(itrial,nCyc(itrial)-1);
        end
    end
    
    rad_targ_miss_all(1,imouse) = nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,MIx),1),2);
    rad_targ_hit_all(1,imouse) = nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,SIx),1),2);
    
    MIxOSIx = find(MIx+SIx);
    for ioff = 1:noff
        ind = intersect(MIxOSIx,find(targetFramesOff == offs(ioff)));
        rad_targ_isi_all(ioff,imouse) = nanmean(nanmean(radius_target_maxnorm(pre_frames-pre_win:pre_frames-1,ind),1),2);
    end
end

%% figures
figure; 
subplot(2,1,1)
scatter(ones(1,size(mouse_mat,1)), rad_targ_hit_all, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
hold on
errorbar(1, mean(rad_targ_hit_all,2), std(rad_targ_hit_all,[],2)./sqrt(size(mouse_mat,1)), 'ok')
scatter(2.*ones(1,size(mouse_mat,1)), rad_targ_miss_all, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
errorbar(2, mean(rad_targ_miss_all,2), std(rad_targ_miss_all,[],2)./sqrt(size(mouse_mat,1)), 'or')
rad_targ_hit_all_avg = mean(rad_targ_hit_all,2);
rad_targ_hit_all_sem = std(rad_targ_hit_all,[],2)./sqrt(size(mouse_mat,1));
rad_targ_miss_all_avg = mean(rad_targ_miss_all,2);
rad_targ_miss_all_sem = std(rad_targ_miss_all,[],2)./sqrt(size(mouse_mat,1));
ylim([0 1])
ylabel('% of Max Radius')
xlabel('Outcome')
xlim([0 3])
[h, p_rad_targ_outcome] = ttest2(rad_targ_hit_all,rad_targ_miss_all);
title(['p = ' num2str(chop(p_rad_targ_outcome,2))])

rad_targ_isi_all_avg = zeros(1,noff);
rad_targ_isi_all_sem = zeros(1,noff);
subplot(2,1,2)
for ioff = 1:noff
    scatter(offs(ioff).*(1000/frame_rate).*ones(1,size(mouse_mat,1)), rad_targ_isi_all(ioff,:), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5])
    hold on
    errorbar(offs(ioff).*(1000/frame_rate), mean(rad_targ_isi_all(ioff,:),2), std(rad_targ_isi_all(ioff,:),[],2)./sqrt(size(mouse_mat,1)), 'ok')
    rad_targ_isi_all_avg(1,ioff) = mean(rad_targ_isi_all(ioff,:),2);
    rad_targ_isi_all_sem(1,ioff) = std(rad_targ_isi_all(ioff,:),[],2)./sqrt(size(mouse_mat,1));
end    
ylim([0 1])
ylabel('% of Max Radius')
xlabel('ISI (ms)')
xlim([0 1000])
p_rad_targ_isi = anova1(rad_targ_isi_all',[],'off');
title(['p = ' num2str(chop(p_rad_targ_isi,2))])

print(fullfile(LG_base, 'Analysis\2P\Adaptation\preTargRadByOutByISISummary.pdf'),'-dpdf','-bestfit')
save(fullfile(LG_base, 'Analysis\2P\Adaptation\eyeTrackingSummary.mat'), 'mouse_mat', 'date_mat', 'run_str_mat', 'p_rad_targ_outcome', 'p_rad_targ_isi', 'rad_targ_hit_all_avg', 'rad_targ_hit_all_sem',  'rad_targ_miss_all_avg', 'rad_targ_miss_all_sem', 'rad_targ_isi_all_avg', 'rad_targ_isi_all_sem');


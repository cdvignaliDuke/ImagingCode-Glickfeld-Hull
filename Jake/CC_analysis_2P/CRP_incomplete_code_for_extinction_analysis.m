


c2_mean_dfof = squeeze(mean(targetAligndFoverF(:,[maskCat==2],:),2));
c1_mean_dfof = squeeze(mean(targetAligndFoverF(:,[maskCat==1],:),2));
step_per_trial = 1/(size(targetAligndFoverF,3)/3+1);
curr_rgb = [0,0,1];
figure;
for trial_num= 1:3:size(targetAligndFoverF,3)
    
    subplot(1,2,1); hold on;
    plot(tt, c2_mean_dfof(:,trial_num), 'Color', curr_rgb);
    subplot(1,2,2); hold on;
    plot(tt, c1_mean_dfof(:,trial_num), 'Color', curr_rgb);
    
    curr_rgb(1) = curr_rgb(1)+step_per_trial;
    curr_rgb(3) = curr_rgb(3)-step_per_trial;
end




c2_mean_spk = squeeze(mean(targetAlign_events(:,[maskCat==2],:),2));
c1_mean_spk = squeeze(mean(targetAlign_events(:,[maskCat==1],:),2));
trials_per_3rd = [1:floor(size(targetAligndFoverF,3)/3)-1];
aa = length(trials_per_3rd);
aa=55;

figure; 
plot_num=1;
for trial_grp=1:3
    subplot(3,2,plot_num);
    shadedErrorBar(tt, nanmean(c1_mean_spk(:,[trials_per_3rd])*33,2), ...
    nanstd(c1_mean_spk(:,[trials_per_3rd])*33,[],2).\sqrt(aa)  );
    plot_num= plot_num+1;
    ylim([-5 10]);
    
    subplot(3,2,plot_num);
    shadedErrorBar(tt, nanmean(c2_mean_spk(:,[trials_per_3rd])*33,2), ...
    nanstd(c2_mean_spk(:,[trials_per_3rd])*33,[],2).\sqrt(aa)  );
    plot_num= plot_num+1;
    ylim([-5 10]);

    trials_per_3rd = trials_per_3rd+length(trials_per_3rd);
end
[peak_resp_all, peak_resp_ind] = max(delta_resp_all,[],2);
peak_resp_reset = peak_resp_ind;
peak_resp_reset(find(peak_resp_ind == 5)) = 3;
peak_resp_reset(find(peak_resp_ind == 6)) = 2;
peak_resp_reset(find(peak_resp_ind == 7)) = 1;
peak_resp_reset(find(peak_resp_ind == 8)) = 0;

[n, bin] = histc(pref_ori_all_diff(:,3),[0 15 75 91]);

figure; 
subplot 231
scatter(pref_ori_all_diff(good_ind_theta,3), pref_ori_all_diff(good_ind_theta,3)-pref_ori_all_diff(good_ind_theta,1))
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-60 60])
title('250 ms ISI')

ind = intersect(find(theta_90_all(:,1)<22.5),good_ind_theta);
length(ind)
subplot 232
scatter(pref_ori_all_diff(ind,3), pref_ori_all_diff(ind,3)-pref_ori_all_diff(ind,1))
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-60 60])
title('Good fits')

subplot 233
for i = 1:length(n)-1
    ind_pref = intersect(find(bin==i),ind);
    errorbarxy(mean(pref_ori_all_diff(ind_pref,3),1),mean(pref_ori_all_diff(ind_pref,3)-pref_ori_all_diff(ind_pref,1),1),std(pref_ori_all_diff(ind_pref,3),[],1)./sqrt(length(ind_pref)), std(pref_ori_all_diff(ind_pref,3)-pref_ori_all_diff(ind_pref,1),[],1)./sqrt(length(ind_pref))) 
    hold on
end
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-15 15])
xlim([0 90])
title('Good fits')

subplot 234
scatter(pref_ori_all_diff(good_ind_theta,3), pref_ori_all_diff(good_ind_theta,3)-pref_ori_all_diff(good_ind_theta,2))
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-60 60])
title('750 ms ISI')

ind2 = intersect(find(theta_90_all(:,2)<22.5),good_ind_theta);
length(ind2)
subplot 235
scatter(pref_ori_all_diff(ind2,3), pref_ori_all_diff(ind2,3)-pref_ori_all_diff(ind2,2))
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-60 60])
title('Good fits')

subplot 236
for i = 1:length(n)-1
    ind_pref = intersect(find(bin==i),ind2);
    errorbarxy(mean(pref_ori_all_diff(ind_pref,3),1),mean(pref_ori_all_diff(ind_pref,3)-pref_ori_all_diff(ind_pref,2),1),std(pref_ori_all_diff(ind_pref,3),[],1)./sqrt(length(ind_pref)),std(pref_ori_all_diff(ind_pref,3)-pref_ori_all_diff(ind_pref,2),[],1)./sqrt(length(ind_pref))) 
    hold on
end
xlabel('Pref ori (deg)')
ylabel('Peak shift (deg)')
ylim([-15 15])
xlim([0 90])
title('Good fits')

print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'adaptDiffShift_scatter.pdf'),'-dpdf','-fillpage')

%%
ind_orthog = intersect(find(~isnan(delta_resp_all(:,1,3))),find(ori_sig_all(:,4)+ori_sig_all(:,8)==2));
figure; 
subplot 221
scatter(delta_resp_all(ind_orthog, 4, 3),delta_resp_all(ind_orthog, 4, 1),[],delta_resp_all(ind_orthog, 8, 3))
xlim([-0.1 1])
ylim([-0.1 1])
refline(1,0)
xlabel('dF/F - 90 deg - Control')
ylabel('dF/F - 90 deg - 250 ms ISI')
[h,p] = ttest(delta_resp_all(ind_orthog, 4, 3),delta_resp_all(ind_orthog, 4, 1));
[r,p] = corrcoef(delta_resp_all(ind_orthog, 4, 3),delta_resp_all(ind_orthog, 4, 1));
title(['p = ' num2str(chop(p,2))])
subplot 222
hist((delta_resp_all(ind_orthog, 4, 3)-delta_resp_all(ind_orthog, 4, 1))./(delta_resp_all(ind_orthog, 4, 3)+delta_resp_all(ind_orthog, 4, 1)))
xlabel('(Control-250)/(Control+250)')
vline(nanmean((delta_resp_all(ind_orthog, 4, 3)-delta_resp_all(ind_orthog, 4, 1))./(delta_resp_all(ind_orthog, 4, 3)+delta_resp_all(ind_orthog, 4, 1)),1))
title(['mean = ' num2str(chop(nanmean((delta_resp_all(ind_orthog, 4, 3)-delta_resp_all(ind_orthog, 4, 1))./(delta_resp_all(ind_orthog, 4, 3)+delta_resp_all(ind_orthog, 4, 1)),1),2))])
subplot 223
scatter(delta_resp_sub_all(ind_orthog, 4, 3),delta_resp_sub_all(ind_orthog, 4, 1),[],delta_resp_sub_all(ind_orthog, 8, 3))
xlim([-0.1 1])
ylim([-0.1 1])
refline(1,0)
xlabel('dF/F - 90 deg - Control')
ylabel('dF/F - 90 deg - 250 ms ISI')
[h,p] = ttest(delta_resp_sub_all(ind_orthog, 4, 3),delta_resp_sub_all(ind_orthog, 4, 1));
title(['Subtracted- p = ' num2str(chop(p,2))])
subplot 224
hist((delta_resp_all(ind_orthog, 4, 3)-delta_resp_sub_all(ind_orthog, 4, 1))./(delta_resp_sub_all(ind_orthog, 4, 3)+delta_resp_sub_all(ind_orthog, 4, 1)))
xlabel('(Control-250)/(Control+250)')
vline(nanmean((delta_resp_sub_all(ind_orthog, 4, 3)-delta_resp_sub_all(ind_orthog, 4, 1))./(delta_resp_sub_all(ind_orthog, 4, 3)+delta_resp_sub_all(ind_orthog, 4, 1)),1))
title(['Subtracted- mean = ' num2str(chop(nanmean((delta_resp_sub_all(ind_orthog, 4, 3)-delta_resp_sub_all(ind_orthog, 4, 1))./(delta_resp_sub_all(ind_orthog, 4, 3)+delta_resp_sub_all(ind_orthog, 4, 1)),1),2))])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', 'linearity_scatter.pdf'),'-dpdf','-fillpage')



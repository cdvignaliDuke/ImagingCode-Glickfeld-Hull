% sort trials by response at end of trial, across both trial types
v_ci = v(:,ci);
a_ci = a(:,ci);
av_ci = av(:,ci);
[~,sort_ind] = sort(mean(av(late_win,ci),1));
sort_ind = flip(sort_ind);

v_sort = v_ci(:,sort_ind);
a_sort = a_ci(:,sort_ind);
av_sort = av_ci(:,sort_ind);

v_early_resp = mean(v_ci(early_win,:),1);
a_early_resp = mean(a_ci(early_win,:),1);
v_late_resp = mean(v_ci(late_win,:),1);
a_late_resp = mean(a_ci(late_win,:),1);
[~,early_ttest] = ttest(v_early_resp,a_early_resp,'tail','both');
[~,late_ttest] = ttest(v_late_resp,a_late_resp,'tail','both');
[~,early_kstest] = kstest2(v_early_resp,a_early_resp,'tail','unequal');
[~,late_kstest] = kstest2(v_late_resp,a_late_resp,'tail','unequal');

%%
setFigParams4Print('landscape');figure; 
suptitle({'repsonse to trial type - task responsive cells'; 'sorted by late resp (to all trials)'})
colormap(brewermap([],'*RdBu'));

%heatmaps
subplot 131
h = imagesc(v_sort');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('cell #')
title('vis trials');
colorbar
caxis(resp_hm_lim)
subplot 132
h = imagesc(a_sort');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('cell #')
title('aud trials');
colorbar
caxis(resp_hm_lim)
subplot 133
h = imagesc(av_sort');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('cell #')
title('all trials');
colorbar
caxis(resp_hm_lim)

print([fnout '_hm_vis_aud'],'-dpdf','-fillpage')

setFigParams4Print('landscape');figure; 
suptitle({'repsonse to trial type - task responsive cells'})
colormap(brewermap([],'*RdBu'));
%time-course
subplot 131
p = patch(early_win_patch(:,1), early_win_patch(:,2), [0.75 0.75 0.75]);
p.EdgeColor = [1 1 1];
p = patch(late_win_patch(:,1), late_win_patch(:,2), [0.75 0.75 0.75]);
p.EdgeColor = [1 1 1];
hold on
h = shadedErrorBar([],mean(v_ci,2),ste(v_ci,2),'-g');
h.edge = [1 1 1];
h.mainLine = 2;
hold on
shadedErrorBar([],mean(a_ci,2),ste(a_ci,2),'k-');
h.edge = [1 1 1];
h.mainLine = 2;
hold on
figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
figYAxis([],'dF/F',resp_y_lim);
figAxForm([]);
vline(bl_fr,'k--');

subplot 132
h = cdfplot(v_early_resp);
h.Color = 'g';
hold on
h = cdfplot(a_early_resp);
h.Color = 'k';
figXAxis([],'time (s)',resp_cdf_lim);
figYAxis([],'dF/F',[0 1]);
figAxForm([]);
title({'early win'; ['test,p = ' num2str(early_ttest)]})

subplot 133
h = cdfplot(v_late_resp);
h.Color = 'g';
hold on
h = cdfplot(a_late_resp);
h.Color = 'k';
figXAxis([],'time (s)',resp_cdf_lim);
figYAxis([],'dF/F',[0 1]);
figAxForm([]);
title({'late win'; ['test,p = ' num2str(late_ttest)]})

print([fnout '_tc_avQuant'],'-dpdf','-fillpage')
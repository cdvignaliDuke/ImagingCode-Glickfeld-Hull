vs = vTC_shuff(tr_start:end,cell_ind);
as = aTC_shuff(tr_start:end,cell_ind);

vs_tc = nanmean(vs(:,ci),2);
vs_err = ste(vs(:,ci),2);
as_tc = nanmean(as(:,ci),2);
as_err = ste(as(:,ci),2);

vsuba_shuff = vs-as;
vsuba_shuff_ci = vsuba_shuff(:,ci);
[~,sort_ind_shuff] = sort(mean(vsuba_shuff_ci(late_win,:),1));
sort_ind_shuff = fliplr(sort_ind_shuff);
vsuba_shuff_sort = vsuba_shuff_ci(:,sort_ind_shuff);
vsuba_shuff_tc = mean(vsuba_shuff(:,ci),2);
vsuba_shuff_err = ste(vsuba_shuff(:,ci),2);

%%
figure; setFigParams4Print('landscape')
suptitle('shuffled')
colormap(brewermap([],'*RdBu'));

subplot 221
h = shadedErrorBar([],vs_tc,vs_err,'-g');
h.edge = [1 1 1];
h.mainLine = 2;
hold on
shadedErrorBar([],as_tc,as_err,'k-');
h.edge = [1 1 1];
h.mainLine = 2;
hold on
figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
figYAxis([],'dF/F',resp_y_lim);
figAxForm([]);
vline(bl_fr,'k--');

subplot 223
h = imagesc(vsuba_shuff_sort');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('all cells')
title({'vis - aud'; 'sorted by late resp'});
colorbar
caxis(sub_hm_lim)

subplot 224
p = patch(early_win_patch(:,1), early_win_patch(:,2), [0.75 0.75 0.75]);
p.EdgeColor = [1 1 1];
p = patch(late_win_patch(:,1), late_win_patch(:,2), [0.75 1 0.75]);
p.EdgeColor = [1 1 1];
hold on
h = shadedErrorBar([],vsuba_shuff_tc,vsuba_shuff_err,'-k');
h.edge = [1 1 1];
h.mainLine = 2;
hold on
h = hline(0,'k--');
vline(bl_fr,'k--');
figXAxis(h.Parent,'time (s)',[1 L],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',sub_y_lim);
figAxForm(h.Parent);
title('vis-aud time-course')

print([fnout '_shuffled'],'-dpdf','-fillpage')

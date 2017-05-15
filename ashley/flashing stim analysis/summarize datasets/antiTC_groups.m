y_axis_lim = [-0.02 0.07];

trL_bin = nbins-3;
vR_long = vR_bin{trL_bin};
aR_long = aR_bin{trL_bin};

vTC_long = vTC_bin{trL_bin};
aTC_long = aTC_bin{trL_bin};



%% responsive cells
figure; setFigParams4Print('landscape');
suptitle(['trial length bin ' num2str(trL_bin)])
subplot(2,2,1)
title('all responsive cells')
hold on
cell_ind = logical(tar_resp | base1_resp | base_sust) & cyc100;
v = vTC_long(tr_start:end,cell_ind);
a = aTC_long(tr_start:end,cell_ind);

% params for linear regression
t_fr = (-bl_fr:size(v,1)-bl_fr-1)';
lr_win = trans_win(end)+1-20:size(v,1);

v_lr = slr(t_fr,nanmean(v,2),lr_win);
a_lr = slr(t_fr,nanmean(a,2),lr_win);

h = plot(t_fr,nanmean(v,2),'g');
h.LineWidth = 1;
hold on
h = plot(t_fr,nanmean(a,2),'k');
h.LineWidth = 1;
hold on
lr = plot(t_fr(lr_win),v_lr,'g--');
lr.LineWidth = 2;
hold on
lr = plot(t_fr(lr_win),a_lr,'k--');
lr.LineWidth = 2;
hold on
vline(0);
figXAxis(h.Parent,'time (s)',[-bl_fr size(v,1)-bl_fr-1],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)

subplot(2,2,2)
title('target resp')
hold on
v = vTC_long(tr_start:end,logical(tar_resp) & cyc100);
a = aTC_long(tr_start:end,logical(tar_resp) & cyc100);
v_lr = slr(t_fr,nanmean(v,2),lr_win);
a_lr = slr(t_fr,nanmean(a,2),lr_win);


h = plot(t_fr,nanmean(v,2),'g');
h.LineWidth = 1;
hold on
h = plot(t_fr,nanmean(a,2),'k');
h.LineWidth = 1;
hold on
lr = plot(t_fr(lr_win),v_lr,'g--');
lr.LineWidth = 2;
hold on
lr = plot(t_fr(lr_win),a_lr,'k--');
lr.LineWidth = 2;
hold on
vline(0);
figXAxis(h.Parent,'time (s)',[-bl_fr size(v,1)-bl_fr-1],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)

subplot(2,2,3)
title('1st base resp')
hold on

v = vTC_long(tr_start:end,logical(base1_resp) & cyc100);
a = aTC_long(tr_start:end,logical(base1_resp) & cyc100);
v_lr = slr(t_fr,nanmean(v,2),lr_win);
a_lr = slr(t_fr,nanmean(a,2),lr_win);

h = plot(t_fr,nanmean(v,2),'g');
h.LineWidth = 1;
hold on
h = plot(t_fr,nanmean(a,2),'k');
h.LineWidth = 1;
hold on
lr = plot(t_fr(lr_win),v_lr,'g--');
lr.LineWidth = 2;
hold on
lr = plot(t_fr(lr_win),a_lr,'k--');
lr.LineWidth = 2;
hold on
vline(0);
figXAxis(h.Parent,'time (s)',[-bl_fr size(v,1)-bl_fr-1],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)

subplot(2,2,4)
title('sustained base resp')
hold on

v = vTC_long(tr_start:end,logical(base_sust) & cyc100);
a = aTC_long(tr_start:end,logical(base_sust) & cyc100);
v_lr = slr(t_fr,nanmean(v,2),lr_win);
a_lr = slr(t_fr,nanmean(a,2),lr_win);

h = plot(t_fr,nanmean(v,2),'g');
h.LineWidth = 1;
hold on
h = plot(t_fr,nanmean(a,2),'k');
h.LineWidth = 1;
hold on
lr = plot(t_fr(lr_win),v_lr,'g--');
lr.LineWidth = 2;
hold on
lr = plot(t_fr(lr_win),a_lr,'k--');
lr.LineWidth = 2;
hold on
vline(0);
figXAxis(h.Parent,'time (s)',[-bl_fr size(v,1)-bl_fr-1],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)


print([fnout '_tc_taskresp'],'-dpdf','-fillpage')
%% auROC groups
figure; setFigParams4Print('landscape');
subplot(2,2,1)
title('auROC-tar2first: decrease')
hold on

v = vTC_long(tr_start:end,logical(auc1_dec & rst1));
a = aTC_long(tr_start:end,logical(auc1_dec & rst1));

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,2,2)
title('auROC-tar2first: increase')
hold on

v = vTC_long(tr_start:end,logical(auc1_inc & rst1));
a = aTC_long(tr_start:end,logical(auc1_inc & rst1));

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,2,3)
title('auROC-tar2last: decrease')
hold on

v = vTC_long(tr_start:end,logical(aucL_dec & rstL));
a = aTC_long(tr_start:end,logical(aucL_dec & rstL));

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,2,4)
title('auROC-tar2last: increase')
hold on

v = vTC_long(tr_start:end,logical(aucL_inc & rstL));
a = aTC_long(tr_start:end,logical(aucL_inc & rstL));

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

print([fnout '_tc_auROC'],'-dpdf','-fillpage')
%% ori pref
figure; setFigParams4Print('landscape');
subplot(2,3,1)
title('0 deg')
hold on

v = vTC_long(tr_start:end,ori_pref == 1);
a = aTC_long(tr_start:end,ori_pref == 1);

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,3,2)
title('45 deg')
hold on

v = vTC_long(tr_start:end,ori_pref == 2);
a = aTC_long(tr_start:end,ori_pref == 2);

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,3,3)
title('90 deg')
hold on

v = vTC_long(tr_start:end,ori_pref == 3);
a = aTC_long(tr_start:end,ori_pref == 3);

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,3,4)
title('135 deg')
hold on

v = vTC_long(tr_start:end,ori_pref == 4);
a = aTC_long(tr_start:end,ori_pref == 4);

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);

subplot(2,3,5)
title('untuned')
hold on

v = vTC_long(tr_start:end,logical(untuned));
a = aTC_long(tr_start:end,logical(untuned));

h = plot(nanmean(v,2),'g');
h.LineWidth = 3;
hold on
h = plot(nanmean(a,2),'k');
h.LineWidth = 3;
hold on
figXAxis(h.Parent,'time (s)',[1 size(v,1)],tr_tick_fr,tr_tick_s);
figYAxis(h.Parent,'dF/F',y_axis_lim);
figAxForm(h.Parent)
hold on
vline(bl_fr);
print([fnout '_tc_tuned'],'-dpdf','-fillpage')

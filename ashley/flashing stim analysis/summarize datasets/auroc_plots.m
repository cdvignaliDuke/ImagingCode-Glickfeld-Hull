y_axis_lim = [-0.04 0.06];
resp_axis_lim = [-0.2 0.2];

figure;
histogram(allAuROC_all_sort,20);hold on;
h = histogram(allAuROC_all_sort(allRst_all_sort),20);
figAxForm(h.Parent);
figXAxis(h.Parent,'auROC',[0 1])
figYAxis(h.Parent,'n cells',[])
legend({'all cells'; 'signif'})
title('all trials')

[auROC_inc_all,auROC_dec_all,auROC_none_all,auROC_grp_all] = auROCsignifGroups(allAuROC_all_sort,allRst_all_sort);
[auROC_inc_short,auROC_dec_short,auROC_none_short,auROC_grp_short] = auROCsignifGroups(allAuROC_short_sort,allRst_short_sort);
[auROC_inc_long,auROC_dec_long,auROC_none_long,auROC_grp_long] = auROCsignifGroups(allAuROC_long_sort,allRst_long_sort);

[hitAuROC_inc_all,hitAuROC_dec_all,hitAuROC_none_all,hitAuROC_grp_all] = auROCsignifGroups(hitAuROC_all_sort,hitRst_all_sort);
[hitAuROC_inc_short,hitAuROC_dec_short,hitAuROC_none_short,hitAuROC_grp_short] = auROCsignifGroups(hitAuROC_short_sort,hitRst_short_sort);
[hitAuROC_inc_long,hitAuROC_dec_long,hitAuROC_none_long,hitAuROC_grp_long] = auROCsignifGroups(hitAuROC_long_sort,hitRst_long_sort);

[missAuROC_inc_all,missAuROC_dec_all,missAuROC_none_all,missAuROC_grp_all] = auROCsignifGroups(missAuROC_all_sort,missRst_all_sort);
[missAuROC_inc_short,missAuROC_dec_short,missAuROC_none_short,missAuROC_grp_short] = auROCsignifGroups(missAuROC_short_sort,missRst_short_sort);
[missAuROC_inc_long,missAuROC_dec_long,missAuROC_none_long,missAuROC_grp_long] = auROCsignifGroups(missAuROC_long_sort,missRst_long_sort);

auroc_color = {'r';'b';'k'};

%% plot timecourses of each auROC group for valid and invalid trials
auroc_plots_allTrL

%% organize data and labelling 
auroc_title = {'inc';'dec';'none'};
trL_title = {'all';'short';'long'};
axes = {x_axis_lim,y_axis_lim,resp_axis_lim};
titles = {auroc_title;trL_title};
tc_ax = {ttMs;trL_ind;tr_tick_s};

%% plot tcs for all, long, and short trials

% all outcomes
d_val = {allVal_all_sort;allVal_long_sort;allVal_short_sort};
d_inv = {allInv_all_sort;allInv_long_sort;allInv_short_sort};
auroc_ind = {auROC_grp_all;auROC_grp_long;auROC_grp_short};
[aurocTC_trL_fig, aurocCDF_trL_fig, aurocScat_trL_fig, aurocHist_trL_fig] = auroc_plots_grpXtrL(d_val,d_inv,auroc_ind,trans_win,titles,axes,tc_ax);
figure(aurocTC_trL_fig)
print([fnout titleStr '_auroc_tc'],'-dpdf','-fillpage')
figure(aurocCDF_trL_fig)
print([fnout titleStr '_auroc_cdf'],'-dpdf','-fillpage')
figure(aurocScat_trL_fig)
print([fnout titleStr '_auroc_scat'],'-dpdf','-fillpage')
figure(aurocHist_trL_fig)
print([fnout titleStr '_auroc_hist'],'-dpdf','-fillpage')

% hits only
d_val = {hitVal_all_sort;hitVal_long_sort;hitVal_short_sort};
d_inv = {hitInv_all_sort;hitInv_long_sort;hitInv_short_sort};
auroc_ind = {hitAuROC_grp_all;hitAuROC_grp_long;hitAuROC_grp_short};
[aurocTC_trL_fig, aurocCDF_trL_fig, aurocScat_trL_fig, aurocHist_trL_fig] = auroc_plots_grpXtrL(d_val,d_inv,auroc_ind,trans_win,titles,axes,tc_ax);
figure(aurocTC_trL_fig)
print([fnout titleStr '_auroc_tc_hits'],'-dpdf','-fillpage')
figure(aurocCDF_trL_fig)
print([fnout titleStr '_auroc_cdf_hits'],'-dpdf','-fillpage')
figure(aurocScat_trL_fig)
print([fnout titleStr '_auroc_scat_hits'],'-dpdf','-fillpage')
figure(aurocHist_trL_fig)
print([fnout titleStr '_auroc_hist_hits'],'-dpdf','-fillpage')

% misses only
d_val = {missVal_all_sort;missVal_long_sort;missVal_short_sort};
d_inv = {missInv_all_sort;missInv_long_sort;missInv_short_sort};
auroc_ind = {missAuROC_grp_all;missAuROC_grp_long;missAuROC_grp_short};
[aurocTC_trL_fig, aurocCDF_trL_fig, aurocScat_trL_fig, aurocHist_trL_fig] = auroc_plots_grpXtrL(d_val,d_inv,auroc_ind,trans_win,titles,axes,tc_ax);
figure(aurocTC_trL_fig)
print([fnout titleStr '_auroc_tc_misses'],'-dpdf','-fillpage')
figure(aurocCDF_trL_fig)
print([fnout titleStr '_auroc_cdf_misses'],'-dpdf','-fillpage')
figure(aurocScat_trL_fig)
print([fnout titleStr '_auroc_scat_misses'],'-dpdf','-fillpage')
figure(aurocHist_trL_fig)
print([fnout titleStr '_auroc_hist_misses'],'-dpdf','-fillpage')
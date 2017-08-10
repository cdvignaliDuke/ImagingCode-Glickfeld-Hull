%% use vis_rep and aud_rep (and rep_str) from plotDepPrevTrial

ind = ci;
vsuba_rep = cellfun(@(x,y) x-y, vis_rep, aud_rep, 'unif', 0);

%% plot change in vis-aud modulation with repeated trials
setFigParams4Print('landscape'); figure;
suptitle({'vis-aud modualtion'; 'red and blue cells mod across all trials'; 'green cells mod across 2 repeats'})
subplot 121
resp1 = mean(vsuba_rep{1}(late_win,ind),1);
resp2 = mean(vsuba_rep{2}(late_win,ind),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = [.75 .75 .75];
h.MarkerEdgeColor = [1,1,1];
hold on
resp1 = mean(vsuba_rep{1}(late_win,ind & mod_bin_ind{1}),1);
resp2 = mean(vsuba_rep{2}(late_win,ind & mod_bin_ind{1}),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'b';
h.MarkerEdgeColor = [1,1,1];
resp1 = mean(vsuba_rep{1}(late_win,ind & mod_bin_ind{2}),1);
resp2 = mean(vsuba_rep{2}(late_win,ind & mod_bin_ind{2}),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'r';
h.MarkerEdgeColor = [1,1,1];
resp1 = mean(vsuba_rep{1}(late_win,ind & H_av_late_2rep'),1);
resp2 = mean(vsuba_rep{2}(late_win,ind & H_av_late_2rep'),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'none';
h.MarkerEdgeColor = 'g';
plot(sub_prev_lim,sub_prev_lim,'k--')
vline(0,'k:')
hline(0,'k:')
figXAxis([],rep_str{1},sub_prev_lim)
figYAxis([],rep_str{2},sub_prev_lim)
figAxForm([])

subplot 122
resp1 = mean(vsuba_rep{2}(late_win,ind),1);
resp2 = mean(vsuba_rep{3}(late_win,ind),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = [.75 .75 .75];
h.MarkerEdgeColor = [1,1,1];
hold on
resp1 = mean(vsuba_rep{2}(late_win,ind & mod_bin_ind{1}),1);
resp2 = mean(vsuba_rep{3}(late_win,ind & mod_bin_ind{1}),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'b';
h.MarkerEdgeColor = [1,1,1];
resp1 = mean(vsuba_rep{2}(late_win,ind & mod_bin_ind{2}),1);
resp2 = mean(vsuba_rep{3}(late_win,ind & mod_bin_ind{2}),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'r';
h.MarkerEdgeColor = [1,1,1];
resp1 = mean(vsuba_rep{2}(late_win,ind & H_av_late_2rep'),1);
resp2 = mean(vsuba_rep{3}(late_win,ind & H_av_late_2rep'),1);
h = scatter(resp1,resp2,50,'o');
h.MarkerFaceColor = 'none';
h.MarkerEdgeColor = 'g';
plot(sub_prev_lim,sub_prev_lim,'k--')
vline(0,'k:')
hline(0,'k:')
figXAxis([],rep_str{2},sub_prev_lim)
figYAxis([],rep_str{3},sub_prev_lim)
figAxForm([])

print([fnout '_prevTr_vsuba'],'-dpdf','-fillpage')

%% non-modualted cells tc
nm = ~(mod_bin_ind{1} | mod_bin_ind{2});
figure;
for irep = 1:3
    vtc = nanmean(vis_rep{irep}(:,ind & nm),2);
    verr = ste(vis_rep{irep}(:,ind & nm),2);
    atc = nanmean(aud_rep{irep}(:,ind & nm),2);
    aerr = ste(aud_rep{irep}(:,ind & nm),2);
    
    subplot(1,3,irep)
    h = shadedErrorBar([],vtc,verr,'g');
    hold on
    h = shadedErrorBar([],atc,aerr,'k');
    figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
    figYAxis([],'dF/F',resp_y_lim);
    figAxForm([]);
    h = vline(bl_fr,'k--');
    title(rep_str{irep})
end
print([fnout '_prevCmpTC_nonMod'],'-dpdf','-fillpage')

%% 2 repeats modulated cells tc
mc = H_av_late_2rep';
figure;
for irep = 1:3
    vtc = nanmean(vis_rep{irep}(:,ind & mc),2);
    verr = ste(vis_rep{irep}(:,ind & mc),2);
    atc = nanmean(aud_rep{irep}(:,ind & mc),2);
    aerr = ste(aud_rep{irep}(:,ind & mc),2);
    
    subplot(1,3,irep)
    h = shadedErrorBar([],vtc,verr,'g');
    hold on
    h = shadedErrorBar([],atc,aerr,'k');
    figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
    figYAxis([],'dF/F',resp_y_lim);
    figAxForm([]);
    h = vline(bl_fr,'k--');
    title(rep_str{irep})
end
print([fnout '_prevCmpTC_2repMod'],'-dpdf','-fillpage')
vsuba_pos_ind = mean(vsuba(late_win,:),1) > 0;
vsuba_neg_ind = mean(vsuba(late_win,:),1) < 0;

%%

slct_axis = [-2 2];

ind = ci;

std_pool = stdPool{analysis_bin}(tr_start:end,ind);

anti_diff = vsuba(:,ind)./std_pool;

v_slctv = v(:,ind)./std_pool;
a_slctv = a(:,ind)./std_pool;

v_slctv_early = mean(v_slctv(early_win,:),1);
v_slctv_late = mean(v_slctv(late_win,:),1);
a_slctv_early = mean(a_slctv(early_win,:),1);
a_slctv_late = mean(a_slctv(late_win,:),1);

thresh_tar = 90;

r_diff = mean(anti_diff(late_win,:),1);

[~,sort_ind] = sort(r_diff);

% bot_ind = sort_ind(1:10);
% top_ind = sort_ind(end-9:end);

fr_500ms = floor(0.5*frRateHz);
%%
vt = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        tars = mouse(imouse).expt(iexp).visTargets;
        [~,tar_ind] = min(abs(tars-thresh_tar));
        d = mouse(imouse).expt(iexp).align(tar_align).av(visual).outcome(hits);
        dt = d.stimResp{tar_ind};
        
        vt = cat(2,vt,mean(dt,3));
    end
end
tr = mean(vt(trans_win,:),1)-mean(vt(pre_win,:),1);
vt_bl = mean(vt(pre_win,ind),1);
vt = vt(tr_start:pre_event_frames+fr_500ms,ind);
vt = bsxfun(@minus, vt,vt_bl);

% [~,sort_ind] = sort(tr(ind));
%%
figure; setFigParams4Print('landscape')
suptitle('selectivity of trial types - task responsive cells')
colormap(brewermap([],'*RdBu'));
subplot 131
h = imagesc(anti_diff(:,sort_ind)');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('all cells')
title('vis - aud; sorted by late');
colorbar
caxis(slct_axis)

subplot 132
hm = v_slctv;
h = imagesc(hm(:,sort_ind)');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('all cells')
title('vis selectivity');
colorbar
caxis(slct_axis)
subplot 133
hm = a_slctv;
h = imagesc(hm(:,sort_ind)');
figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(h.Parent);
ylabel('all cells')
title('aud selectivity');
colorbar
caxis(slct_axis)

print([fnout '_selectivity_hm'],'-dpdf','-fillpage')

figure; setFigParams4Print('landscape')
suptitle('')
colormap(brewermap([],'*RdBu'));
subplot 121
resp_lim = [-2.5 2.5];
plot(resp_lim,resp_lim,'b-')
hold on
h = scatter(v_slctv_early,a_slctv_early,100,'k.');
figXAxis(h.Parent,'vis selectivity',resp_lim);
figYAxis(h.Parent,'aud selectivity',resp_lim);
figAxForm(h.Parent);
hold on
vline(0,'k--');
hline(0,'k--');
title('early resp selectivity');
subplot 122
resp_lim = [-2.5 2.5];
plot(resp_lim,resp_lim,'b-')
hold on
h = scatter(v_slctv_late,a_slctv_late,100,'k.');
figXAxis(h.Parent,'vis selectivity',resp_lim);
figYAxis(h.Parent,'aud selectivity',resp_lim);
figAxForm(h.Parent);
hold on
vline(0,'k--');
hline(0,'k--');
title('late resp selectivity');

print([fnout '_selectivity_scat'],'-dpdf','-fillpage')
% %% plot v and a tc for vis groups (pos and neg subtraction cells)
% auc_ind_name = {'task responsive cells'; 'auROC tar:1st base > 0.5'; 'auROC tar:1stbase < 0.5'};
% auc_ind_mat = {logical(1:nc);auc1_inc ; auc1_dec};
% 
% for ifig = 1:length(auc_ind_mat);
% loop_ind = ci & auc_ind_mat{ifig};
% 
% y_lim = [-0.015 0.03];
% %increase on vis trials
% ind = loop_ind & vsuba_pos_ind;
% 
% clear leg
% figure;
% setFigParams4Print('landscape')
% colormap(brewermap([],'*RdBu'));
% suptitle(auc_ind_name{ifig})
% subplot 221
% v_mtc = mean(v(:,ind),2);
% a_mtc = mean(a(:,ind),2);
% v_err = ste(v(:,ind),2);
% a_err = ste(a(:,ind),2);
% 
% vsuba_pos_mtc = mean(anti_diff(:,ind),2);
% vsuba_pos_err = ste(anti_diff(:,ind),2);
% 
% h = shadedErrorBar([],v_mtc,v_err,'g');
% h.mainLine.LineWidth = 2;
% leg(1) = h.mainLine;
% hold on
% h = shadedErrorBar([],a_mtc,a_err,'k');
% h.mainLine.LineWidth = 2;
% leg(2) = h.mainLine;
% h = vline(bl_fr,'k--');
% figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm(h.Parent);
% figYAxis(h.Parent,'dF/F',y_lim)
% legend(leg,{'vis','aud'},'location','southwest')
% title('cells: v - a > 0')
% 
% %decrease on vis trials
% ind = loop_ind & vsuba_neg_ind;
% 
% subplot 222
% v_mtc = mean(v(:,ind),2);
% a_mtc = mean(a(:,ind),2);
% v_err = ste(v(:,ind),2);
% a_err = ste(a(:,ind),2);
% 
% vsuba_neg_mtc = mean(anti_diff(:,ind),2);
% vsuba_neg_err = ste(anti_diff(:,ind),2);
% 
% h = shadedErrorBar([],v_mtc,v_err,'g');
% h.mainLine.LineWidth = 2;
% leg(1) = h.mainLine;
% hold on
% h = shadedErrorBar([],a_mtc,a_err,'k');
% h.mainLine.LineWidth = 2;
% leg(2) = h.mainLine;
% h = vline(bl_fr,'k--');
% figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm(h.Parent);
% figYAxis(h.Parent,'dF/F',y_lim)
% legend(leg,{'vis','aud'},'location','southwest')
% title('cells: v - a < 0')
% 
% % plot subtractions of each group
% subplot 223
% h = shadedErrorBar([],vsuba_pos_mtc,vsuba_pos_err,'r');
% h.mainLine.LineWidth = 2;
% leg(1) = h.mainLine;
% hold on
% h = shadedErrorBar([],vsuba_neg_mtc,vsuba_neg_err,'b');
% h.mainLine.LineWidth = 2;
% leg(2) = h.mainLine;
% h = vline(bl_fr,'k--');
% figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm(h.Parent);
% figYAxis(h.Parent,'dF/F',[-0.015 0.015])
% legend(leg,{'v-a > 0','v-a < 0'},'location','northwest')
% title('v - a subtraction across groups')
% 
% %plot heatmap of subtraction
% subplot 224
% h = imagesc(flipud(anti_diff(:,sort_ind)'));
% figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm(h.Parent);
% ylabel('all cells')
% title('vis - aud; sorted by late');
% colorbar
% caxis([-0.1 0.1])
% 
% end
% 
% %% plot tcs and subtraction for auROC groups
% 
% y_lim = [-0.015 0.03];
% 
% clear leg
% figure;
% setFigParams4Print('landscape')
% colormap(brewermap([],'*RdBu'));
% suptitle('')
% 
% vsuba_auc_mtc = cell(1,2);
% vsuba_auc_err = cell(1,2);
% for iplot = 1:2
% ind = ci & auc_ind_mat{iplot+1};
% 
% subplot(2,2,iplot)
% v_mtc = mean(v(:,ind),2);
% a_mtc = mean(a(:,ind),2);
% v_err = ste(v(:,ind),2);
% a_err = ste(a(:,ind),2);
% 
% vsuba_auc_mtc{iplot} = mean(anti_diff(:,ind),2);
% vsuba_auc_err{iplot} = ste(anti_diff(:,ind),2);
% 
% h = shadedErrorBar([],v_mtc,v_err,'g');
% h.mainLine.LineWidth = 2;
% leg(1) = h.mainLine;
% hold on
% h = shadedErrorBar([],a_mtc,a_err,'k');
% h.mainLine.LineWidth = 2;
% leg(2) = h.mainLine;
% figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm([]);
% figYAxis([],'dF/F',y_lim)
% h = vline(bl_fr,'k--');
% legend(leg,{'vis','aud'},'location','southwest')
% title(auc_ind_name{iplot+1})
% end
% 
% subplot 223
% h = shadedErrorBar([],vsuba_auc_mtc{1},vsuba_auc_err{1},'r');
% h.mainLine.LineWidth = 2;
% leg(1) = h.mainLine;
% hold on
% h = shadedErrorBar([],vsuba_auc_mtc{2},vsuba_auc_err{2},'b');
% h.mainLine.LineWidth = 2;
% leg(2) = h.mainLine;
% h = vline(bl_fr,'k--');
% figXAxis(h.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
% figAxForm(h.Parent);
% figYAxis(h.Parent,'dF/F',[-0.015 0.015])
% legend(leg,auc_ind_name(2:3),'location','northwest')
% title('v - a subtraction across groups')
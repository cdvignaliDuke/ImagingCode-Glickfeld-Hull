
%% all trial lengths and types
aurocFig = figure;
suptitle('all trial lengths')
d = eval(dn_all_val{3});
subplot(1,2,1)
for itc = 1:3

   cell_ind = auROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('all valid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);
d = eval(dn_all_inv{3});
subplot(1,2,2)
for itc = 1:3

   cell_ind = auROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('all invalid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);
print([fnout titleStr '_auroc_val'],'-dpdf','-fillpage')
%% all trial lengths, hits only
aurocFig = figure;
suptitle('all trial lengths')
d = eval(dn_all_val{1});
subplot(1,2,1)
for itc = 1:3

   cell_ind = hitAuROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('hit valid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);

d = eval(dn_all_inv{1});
subplot(1,2,2)
for itc = 1:3

   cell_ind = hitAuROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('hit invalid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);

%% all trial lengths, misses only
aurocFig = figure;
suptitle('all trial lengths')
d = eval(dn_all_val{2});
subplot(1,2,1)
for itc = 1:3

   cell_ind = missAuROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('miss valid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);

d = eval(dn_all_inv{2});
subplot(1,2,2)
for itc = 1:3

   cell_ind = missAuROC_grp_all{itc};
    tc = mean(d(cell_ind,trL_ind),1);
    tc_err = ste(d(cell_ind,trL_ind),1);
    grptc(itc) = shadedErrorBar(ttMs,tc,tc_err,'k');        
    leg = grptc(itc).mainLine;
    leg.Color = auroc_color{itc};
    leg.LineWidth = 3;
    hold on            
end
title('miss invalid')
figXAxis(grptc(itc),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
figYAxis(grptc(itc),'dF/F',y_axis_lim);
hold on
vl = vline(0,'k--');
figAxForm(vl.Parent);
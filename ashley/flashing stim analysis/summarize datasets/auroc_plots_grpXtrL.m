function [aurocTC_trL_fig, aurocCDF_trL_fig, aurocScat_trL_fig, aurocHist_trL_fig] = auroc_plots_grpXtrL(d_val,d_inv,auroc_ind,trans_win,titles,axes,tc_ax)
x_axis_lim = axes{1};
y_axis_lim = axes{2};
resp_axis_lim = axes{3};

auroc_title = titles{1};
trL_title = titles{2};

ttMs = tc_ax{1};
trL_ind = tc_ax{2};
tr_tick_s = tc_ax{3};

aurocTC_trL_fig = figure;setFigParams4Print('landscape')
aurocCDF_trL_fig = figure;setFigParams4Print('landscape')
aurocScat_trL_fig = figure;setFigParams4Print('landscape')
aurocHist_trL_fig = figure;setFigParams4Print('landscape')

for iplot = 1:3
    p1 = 1 + ((iplot-1)*3);
%     p2 = 2 + ((iplot-1)*3);
%     p3 = 3 + ((iplot-1)*3);
    ind = auroc_ind{iplot};
    val_temp = d_val{iplot};
    inv_temp = d_inv{iplot};
    resp_val_temp = mean(val_temp(:,trans_win),2);
    resp_inv_temp = mean(inv_temp(:,trans_win),2);
    for iind = 1:3
        cell_ind = ind{iind}; % inc, dec, or none
%         time-courses
        figure(aurocTC_trL_fig);
       sp = subplot(3,3,p1+(iind-1));
        tc = mean(val_temp(cell_ind,trL_ind),1);
        tc_err = ste(val_temp(cell_ind,trL_ind),1);
        grptc(1) = shadedErrorBar(ttMs,tc,tc_err,'k');    
        leg = grptc(1).mainLine;
        leg.LineWidth = 3;
        hold on  
        tc = mean(inv_temp(cell_ind,trL_ind),1);
        tc_err = ste(inv_temp(cell_ind,trL_ind),1);
        grptc(2) = shadedErrorBar(ttMs,tc,tc_err,'c');    
        leg = grptc(2).mainLine;
        leg.LineWidth = 3;
        hold on  
        title(auroc_title{iind})
        figXAxis(grptc(2),'time (s)',x_axis_lim,tr_tick_s,tr_tick_s);
        figYAxis(grptc(2),{trL_title{iplot};'dF/F'},y_axis_lim);
        hold on
        vl = vline(0,'k--');
        figAxForm(vl.Parent);
        
%         cdfs
        figure(aurocCDF_trL_fig)
        sp = subplot(3,3,p1+(iind-1));
        l = cdfplot(resp_val_temp(cell_ind));
        l.Color = 'k';
        l.LineWidth = 1;
        hold on
        l = cdfplot(resp_inv_temp(cell_ind));
        l.Color = 'c';
        l.LineWidth = 1;
        title(auroc_title{iind})
        figXAxis(l,'dF/F',resp_axis_lim);
        figYAxis(l,{trL_title{iplot};'% cells'},[0 1]);
        hold on
        vl = vline(0,'k--');
        figAxForm(vl.Parent);      
        
%         scatters
        figure(aurocScat_trL_fig)
        sp = subplot(3,3,p1+(iind-1));
        s = scatter(resp_val_temp(cell_ind),resp_inv_temp(cell_ind),100,'k.');
        title(auroc_title{iind})
        hold on
        vl = plot(-20:1:20,-20:1:20,'k--');
        figXAxis(s,'val dF/F',resp_axis_lim);
        figYAxis(s,{trL_title{iplot};'inv dF/F'},resp_axis_lim);
        figAxForm(vl.Parent);   

%         histogram of b from scatter (y = mx+b, m=1, b = y-x)
        figure(aurocHist_trL_fig)
        b = (resp_val_temp(cell_ind)-resp_inv_temp(cell_ind));
        b_edges = linspace(-max(abs(b)), max(abs(b)),10);
        sp = subplot(3,3,p1+(iind-1));
        h = histogram(b,b_edges);
        title(auroc_title{iind})
        hold on
        figXAxis(h,'(val-inv)',[b_edges(1) b_edges(end)]);
        figYAxis(h,{trL_title{iplot};'n cells'},[0 30]);
        vl = vline(mean(b),'r-');
        hold on
        vl = vline(0,'k--');
        figAxForm(vl.Parent); 
    end
end
end
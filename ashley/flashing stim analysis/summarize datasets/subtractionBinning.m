%% percentile bins
% get shuffed tc
vs = vTC_shuff(tr_start:end,cell_ind);
as = aTC_shuff(tr_start:end,cell_ind);
vssubas = vs-as;
late_sub_shuff = mean(vssubas(late_win,:));

% get indexes for percentile bins
[~,sort_ind_late] = sort(late_sub);
[~,sort_ind_late_shuff] = sort(late_sub_shuff);
nbins = 4;
nc = length(late_sub);
bin_size = floor(nc/nbins);

cell_bin_ind = cell(1,nbins);
cell_bin_shuff_ind = cell(1,nbins);
offset = 0;
for i = 1:nbins
    if i < nbins
        cell_bin_ind{i} = sort_ind_late(offset+1:offset+bin_size);
        cell_bin_shuff_ind{i} = sort_ind_late_shuff(offset+1:offset+bin_size);
        offset = offset+bin_size;
    else
        
        cell_bin_ind{i} = sort_ind_late(offset+1:end);
        cell_bin_shuff_ind{i} = sort_ind_late_shuff(offset+1:end);
    end
end
%% binned mean time-course, all task-responsive cells

vs_bin = cell(1,nbins);
as_bin = cell(1,nbins);

v_bin = cell(1,nbins);
a_bin = cell(1,nbins);

for i = 1:nbins
   v_bin{i} = v(:,intersect(cell_bin_ind{i},find(ci)));
   a_bin{i} = a(:,intersect(cell_bin_ind{i},find(ci))); 
   vs_bin{i} = vs(:,intersect(cell_bin_shuff_ind{i},find(ci)));
   as_bin{i} = as(:,intersect(cell_bin_shuff_ind{i},find(ci))); 
end

v_bin_mean = cellfun(@(x) mean(x,2),v_bin,'unif',0);
a_bin_mean = cellfun(@(x) mean(x,2),a_bin,'unif',0);
vs_bin_mean = cellfun(@(x) mean(x,2),vs_bin,'unif',0);
as_bin_mean = cellfun(@(x) mean(x,2),as_bin,'unif',0);

% colors = brewermap(nbins+1,'*Greens');
colors = brewermap(nbins+1,'*Greys');

figure;
suptitle('percentile bins')
for i = 1:nbins
   subplot(2,2,1)
   h = plot(v_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       y_axis_lim = h.Parent.YLim;
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',[]);
       figAxForm([])
       title('vis trials')
       leg = legend(num2str([1:nbins]'),'location','northeastoutside');
       title(leg,'bin #')
   end
   subplot(2,2,2)
   h = plot(a_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',[y_axis_lim]);
       figAxForm([])
       title('aud trials')
       leg = legend(num2str([1:nbins]'),'location','northeastoutside');
       title(leg,'bin #')
   end
   subplot(2,2,3)
   h = plot(vs_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       y_axis_lim = h.Parent.YLim;
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',y_axis_lim);
       figAxForm([])
       title('vis trials - shuffled')
       leg = legend(num2str([1:nbins]'),'location','northeastoutside');
       title(leg,'bin #')
   end
   subplot(2,2,4)
   h = plot(as_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',[y_axis_lim]);
       figAxForm([])
       title('aud trials - shuffed')
       leg = legend(num2str([1:nbins]'),'location','northeastoutside');
       title(leg,'bin #')
   end
end

print([fnout '_subBins_tc'],'-dpdf','-fillpage')

%% percentile bins
nbins = 2;

cell_bin_ind = cell(1,nbins);

cell_bin_ind{1} = late_sub < 0;
cell_bin_ind{2} = late_sub > 0;


%% binned mean time-course, significantly modulated cells
skipnans = ~isnan(H_av_late);
v_bin = cell(1,nbins);
a_bin = cell(1,nbins);
for i = 1:nbins
   v_bin{i} = v(:,cell_bin_ind{i}(skipnans) & H_av_late(skipnans) & ci(skipnans));
   a_bin{i} = a(:,cell_bin_ind{i}(skipnans) & H_av_late(skipnans) & ci(skipnans)); 
end

v_bin_mean = cellfun(@(x) mean(x,2),v_bin,'unif',0);
a_bin_mean = cellfun(@(x) mean(x,2),a_bin,'unif',0);

% colors = brewermap(nbins+1,'*Greens');
colors = brewermap(nbins+1,'*Greys');
nc_bin = cellfun(@(x) size(x,2),v_bin);

figure;
suptitle('significantly modulated cells')
for i = 1:nbins
   subplot(1,2,1)
   h = plot(v_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   leg_str{i} = sprintf('%s cells', num2str(nc_bin(i)));
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',y_axis_lim);
       figAxForm([])
       title('vis trials')
       leg = legend(leg_str,'location','northeastoutside');
       title(leg,'bin #')
   y_axis_lim = h.Parent.YLim;
   end
   subplot(1,2,2)
   h = plot(a_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',y_axis_lim);
       figAxForm([])
       title('aud trials')
       leg = legend(leg_str,'location','northeastoutside');
       title(leg,'bin #')
   end
end

print([fnout '_modBins_tc'],'-dpdf','-fillpage')

%% histogram of cell types in each bin
figure;
suptitle('task-responsive cells')
colormap gray
% orientation tuning
nOriCells_bin = zeros(nbins,5);
for ibin = 1:nbins
    for iori = 1:5
        if iori == 5
            ind = logical(untuned) & ci;
        else
            ind = ori_pref == iori & ci;
        end
        nOriCells_bin(ibin,iori) = sum(ind(cell_bin_ind{ibin}));
    end
end
pctOriCells_bin = bsxfun(@rdivide, nOriCells_bin,cellfun(@(x) sum(x & ci),cell_bin_ind)');

subplot 121
bar(1:nbins,pctOriCells_bin,'grouped')
title('orientation tuning')
leg = legend({'0';'45';'90';'135';'untuned'},'location','northeastoutside');
title(leg,'ori tuning')
figXAxis([], 'late subtraction bin',[],1:nbins,{'vis<aud';'vis>aud'})
figYAxis([],'% cells in bin',[])
figAxForm([])

if ~strcmp('_audControl',ds)
    % auROC baseline:target
    nAucCells_bin = zeros(nbins,2);
    for ibin = 1:nbins
        ind = auc1_dec & rst1_dir & ci;
        nAucCells_bin(ibin,1) = sum(ind(cell_bin_ind{ibin}));
        ind = auc1_inc & rst1_dir & ci;
        nAucCells_bin(ibin,2) = sum(ind(cell_bin_ind{ibin}));
    end
    pctAucCells_bin = bsxfun(@rdivide, nAucCells_bin,cellfun(@(x) sum(x & ci),cell_bin_ind)');

    subplot 122
    bar(1:nbins,pctAucCells_bin,'grouped')
    title('auROC: 1st BL stim to Target stim')
    leg = legend({'bl>tar';'bl<tar'},'location','northeastoutside');
    title(leg,'auROC result')
    figXAxis([], 'late subtraction bin',[],1:nbins,{'vis<aud';'vis>aud'})
    figYAxis([],'% cells in bin',[])
    figAxForm([])

    print([fnout '_subBins_taskResp'],'-dpdf','-fillpage')

    figure;
    suptitle('late window modulated')
    colormap gray
    % orientation tuning
    nOriCells_bin = zeros(nbins,5);
    for ibin = 1:nbins
        for iori = 1:5
            if iori == 5
                ind = logical(untuned) & H_av_late;
            else
                ind = ori_pref == iori & H_av_late;
            end
            nOriCells_bin(ibin,iori) = sum(ind(cell_bin_ind{ibin}));
        end
    end
    pctOriCells_bin = bsxfun(@rdivide, nOriCells_bin,cellfun(@(x) sum(x & H_av_late),cell_bin_ind)');

    subplot 121
    bar(1:nbins,pctOriCells_bin,'grouped')
    title('orientation tuning')
    leg = legend({'0';'45';'90';'135';'untuned'},'location','northeastoutside');
    title(leg,'ori tuning')
    figXAxis([], 'late subtraction bin',[],1:nbins,{'vis<aud';'vis>aud'})
    figYAxis([],'% cells in bin',[])
    figAxForm([])

    % auROC baseline:target
    nAucCells_bin = zeros(nbins,2);
    for ibin = 1:nbins
        ind = auc1_dec & rst1_dir & H_av_late;
        nAucCells_bin(ibin,1) = sum(ind(cell_bin_ind{ibin}));
        ind = auc1_inc & rst1_dir & H_av_late;
        nAucCells_bin(ibin,2) = sum(ind(cell_bin_ind{ibin}));
    end
    pctAucCells_bin = bsxfun(@rdivide, nAucCells_bin,cellfun(@(x) sum(x & H_av_late),cell_bin_ind)');

    subplot 122
    bar(1:nbins,pctAucCells_bin,'grouped')
    title('auROC: 1st BL stim to Target stim')
    leg = legend({'bl>tar';'bl<tar'},'location','northeastoutside');
    title(leg,'auROC result')
    figXAxis([], 'late subtraction bin',[],1:nbins,{'vis<aud';'vis>aud'})
    figYAxis([],'% cells in bin',[])
    figAxForm([])

    print([fnout '_subBins_vsubaMod'],'-dpdf','-fillpage')
end
%% highlight modulated cells in vis vs aud response scatter plot
v_late = mean(v(late_win,ci),1);
a_late = mean(a(late_win,ci),1);
resp_lim = [-0.15 0.5];
if strcmp('_audControl',ds)
    y_axis_lim = [-0.04 0.065];
else
    y_axis_lim = [-0.025 0.04];
end
leg = [];
setFigParams4Print('landscape'); figure;
suptitle('task-responsive cells')

subplot 131
h = scatter(v_late,a_late,'ko');
h.MarkerFaceColor = [.75 .75 .75];
h.MarkerEdgeColor = [1 1 1];
hold on
ind = logical(H_av_late(ci)) & cell_bin_ind{1}(ci);
h = scatter(v_late(ind),a_late(ind),'ko');
h.MarkerFaceColor = 'b';
h.MarkerEdgeColor = [1 1 1];
leg(1) = h;
ind = logical(H_av_late(ci)) & cell_bin_ind{2}(ci);
h = scatter(v_late(ind),a_late(ind),'ko');
h.MarkerFaceColor = 'r';
h.MarkerEdgeColor = [1 1 1];
leg(2) = h;
plot(resp_lim,resp_lim,'k--')

figXAxis([],'visual',resp_lim);
figYAxis([],'auditory',resp_lim);
figAxForm([])
title('late window response')
legend(leg,{'vis<aud';'vis>aud'},'location','northwest')


for i = 1:nbins
   subplot 132
   h = plot(v_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   leg_str{i} = sprintf('%s cells', num2str(nc_bin(i)));
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',y_axis_lim);
       figAxForm([])
       title('vis trials')
       legend(leg_str,'location','northwest');
       y_axis_lim = h.Parent.YLim;
   end
   subplot 133
   h = plot(a_bin_mean{i});
   h.Color = colors(i,:);
   h.LineWidth = 2;
   hold on
   if i == nbins
       figXAxis([],'time (s)',[],tr_tick_fr,tr_tick_s);
       figYAxis([],'dF/F',y_axis_lim);
       figAxForm([])
       title('aud trials')
   end
end

print([fnout '_scat_modCells'],'-dpdf','-fillpage')

%% dependence on previous trial
clear leg
fig_ttl_str = {'visLTaud';'visGTaud'};
% y_lim = [-0.005 0.02];
% y_lim = [-0.025 0.1];

skipnans = ~isnan(H_av_late);
mod_bin_ind = cell(1,nbins);
for i = 1:nbins
   mod_bin_ind{i} = cell_bin_ind{i}(skipnans) & H_av_late(skipnans) & ci(skipnans);
end

for ibin = 1:nbins
    setFigParams4Print('landscape');figure;
    suptitle(['mod cells, ' fig_ttl_str{ibin}])
    ind = mod_bin_ind{ibin};

    tc1 = mean(v(:,ind),2);
    err1 = ste(v(:,ind),2);
    tc2 = mean(a(:,ind),2);
    err2 = ste(a(:,ind),2);

    subplot 131
    h = shadedErrorBar([],tc1,err1,'g');
    leg(1) = h.mainLine;
    hold on
    h= shadedErrorBar([],tc2,err2,'k');
    leg(2) = h.mainLine;
    figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
    figYAxis([],'dF/F',y_axis_lim);
    figAxForm([]);
    h = vline(bl_fr,'k--');
    legend(leg,{'vis','aud'},'location','northwest')
    title('all trials')

    tc1 = nanmean(vp{1}(:,ind),2);
    err1 = ste(vp{1}(:,ind),2);
    tc2 = nanmean(ap{2}(:,ind),2);
    err2 = ste(ap{2}(:,ind),2);

    subplot 132
    h = shadedErrorBar([],tc1,err1,'g');
    leg(1) = h.mainLine;
    hold on
    h= shadedErrorBar([],tc2,err2,'k');
    leg(2) = h.mainLine;
    figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
    figYAxis([],'dF/F',y_axis_lim);
    figAxForm([]);
    h = vline(bl_fr,'k--');
    legend(leg,{'vis->vis','aud->aud'},'location','northwest')
    title('2 repeats')

    tc1 = nanmean(vp{3}(:,ind),2);
    err1 = ste(vp{3}(:,ind),2);
    tc2 = nanmean(ap{4}(:,ind),2);
    err2 = ste(ap{4}(:,ind),2);

    subplot 133
    h = shadedErrorBar([],tc1,err1,'g');
    leg(1) = h.mainLine;
    hold on
    h = shadedErrorBar([],tc2,err2,'k');
    leg(2) = h.mainLine;
    figXAxis([],'time (s)',[1 L],tr_tick_fr,tr_tick_s);
    figYAxis([],'dF/F',y_axis_lim);
    figAxForm([]);
    h = vline(bl_fr,'k--');
    legend(leg,{'vis-vis-vis';'aud-aud-aud'},'location','northwest')
    title('3 repeats')
    print([fnout '_prevTrCmpTC_' fig_ttl_str{ibin}],'-dpdf','-fillpage')
end
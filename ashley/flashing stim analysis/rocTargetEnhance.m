load([fnout '_modCells'])

trType = 3;

set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
% set(0,'DefaultaxesFontSize', 16)

auc_visTarResp = [];            
auc_visTarResp_short = [];
auc_visTarResp_long = [];
rs_p_visTarResp = []; 
rs_visTarResp = []; 
auc_visTarRespMaxDir = [];
expColor = [];
i=0;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            %auROC across all trials
            a = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value;           
            as = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(1).value;            
            al = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(2).value;
            u = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst_p;
            ut = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst; 
            
            % max auROC across by each target
            b = max(cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value_dirs),[],2);
            
            if size(a,1) > size(a,2)
                a = a';
            end
            auc_visTarResp = cat(2,auc_visTarResp,a);            
            auc_visTarResp_short = cat(1,auc_visTarResp_short,as);
            auc_visTarResp_long = cat(2,auc_visTarResp_long,al);
            rs_visTarResp = cat(2,rs_visTarResp,ut); 
            rs_p_visTarResp = cat(2,rs_p_visTarResp,u); 
            auc_visTarRespMaxDir = cat(1,auc_visTarRespMaxDir,b);
            expColor = cat(1,expColor,ones(length(a),1)+i);
            i = i+1;
        end
    end
end

c = length(auc_visTarResp);
nbins = 4;
val_sort = sort(auc_visTarResp);
percentile_ind = [val_sort(1:floor(c/nbins):c) val_sort(end)];

figure;
clear h
h{1}=cdfplot(auc_visTarResp);
h{1}.Color = 'k';
h{1}.LineWidth = 2;
hold on
vline(percentile_ind(2:end-1),'r--')
hold on
% h{2}=cdfplot(auc_visTarResp);
% h{2}.Color = [0.5 0.5 0.5];
% h{2}.LineWidth = 2;

% h{2}=cdfplot(auc_visTarResp_short);
% h{2}.Color = [0.75 0.75 0.75];
% h{2}.LineWidth = 2;
% 
% h{3}=cdfplot(auc_visTarResp_long);
% h{3}.Color = [0.25 0.25 0.25];
% h{3}.LineWidth = 2;
xlim([0 1])
ylim([0 1])
xlabel('area under ROC curve')
ylabel('fraction of cells')
axis square
title('comp. resp. to vis target and resp. to last vis stim')
% legend({'all';'short';'long'})

print([fnout '_aucCDF'],'-dpdf')

figure;
colormap(jet)
scatter(auc_visTarResp,rs_p_visTarResp,50,expColor);
hold on
hline(0.05,'k--')
xlabel('auc')
ylabel('p value - rank sum test')
colorbar
clim([1 8])

figure;
colormap(brewermap([],'*RdGy'))
clear h
h = scatter(auc_visTarResp,auc_visTarRespMaxDir, 100,rs_visTarResp,'.');
hold on
plot([0:0.1:1],[0:0.1:1],'k--')
xlim([0 1]);
ylim([0 1]);
clim([-0.1 1.1])
colorbar
xlabel('all targets')
ylabel('max across ea. target')
title('auROC')
axis square
% colorbar for ea exp
figure;
clear h
h = scatter(auc_visTarResp,auc_visTarRespMaxDir, 100,expColor,'.');
colormap(jet)
% h = scatter(auc_visTarResp,auc_visTarRespMaxDir, 100,rs_visTarResp,'.');
hold on
plot([0:0.1:1],[0:0.1:1],'k--')
xlim([0 1]);
ylim([0 1]);
% clim([-0.1 1.1])
colorbar
xlabel('all targets')
ylabel('max across ea. target')
title('auROC')
axis square

%% roc values for each target direction

allDirs = [];
c_allDirs = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        allDirs = [allDirs unique(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).dirs)];
        if mouse(imouse).expt(iexp).info.isCatch
        c_allDirs = [c_allDirs unique(msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).dirs)];
        end
    end
end
allDirs = unique(allDirs);
c_allDirs = unique(c_allDirs);

auc_allDirs = cell(1,length(allDirs));
rs_allDirs = [];
rs_tarAbove = [];
rs_tarBelow = [];
resp(hits).all = cell(1,length(allDirs));
resp(fas).all = cell(1,length(c_allDirs));
resp(hits).enhance = cell(1,length(allDirs));
resp(fas).enhance = cell(1,length(c_allDirs));
resp(hits).tuning = cell(3,length(allDirs));
nc = 0;
ncc = 0;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        c = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs')));
        A = cellfun(@(x) find(x > 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs,'unif',false);
        B = cellfun(@(x) find(x < 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs,'unif',false);
        ca = intersect(c,unique(cat(1,cell2mat(A'))));
        cb = intersect(c,unique(cat(1,cell2mat(B'))));
        ind = find(ismember(allDirs,msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).dirs));
%             r = cellfun(@(x,y) x(:,y,:),mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp(2:end),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs,'unif',false);
        for idir = 1:length(ind)
            a = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs{idir};
            auc_allDirs{ind(idir)} = cat(1, auc_allDirs{ind(idir)},a);            
            
            resp(hits).all{ind(idir)} = cat(2, resp(hits).all{ind(idir)}, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp{idir+1},3));
%             resp(hits).enhance{ind(idir)} = cat(2, resp(hits).enhance{ind(idir)}, mean(r{idir},3));
          end
        for ibin = 1:3
            if ibin == 1
                rsDirInd = find(allDirs(ind) < 25);
            elseif ibin == 2                
                rsDirInd = find(allDirs(ind) > 25 & allDirs(ind) < 60);
            elseif ibin == 3
                rsDirInd = find(allDirs(ind) > 60);
            end
            exInd = find(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp,'unif',false)));
            exInd = exInd(2:end);
            rsInd = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs(rsDirInd)')));
            resp(hits).tuning(ibin,ind) = cellfun(@(x,y) cat(2,x,y), resp(hits).tuning(ibin,ind), cellfun(@(z) mean(z(:,rsInd,:),3),mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp(exInd),'unif',false),'unif',false);
        end
        if mouse(imouse).expt(iexp).info.isCatch
        cind = find(ismember(c_allDirs,msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).dirs));
%         trInd = find(~isempty(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp));
%         r = cellfun(@(x,y) x(:,y,:),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp(trInd),msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).rst_dirs,'unif',false);
        for idir = 1:length(cind)
            resp(fas).all{cind(idir)} = cat(2, resp(fas).all{cind(idir)}, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir},3));
%             resp(fas).enhance{cind(idir)} = cat(2, resp(fas).enhance{cind(idir)}, mean(r{idir},3));
        end
        
        end
            %cells index for signif auc for at least 1 target
            
            rs_allDirs = cat(1, rs_allDirs,c+nc);
            rs_tarAbove = cat(1,rs_tarAbove,ca+nc);
            rs_tarBelow = cat(1,rs_tarBelow,cb+nc);
            nc = nc+length(a);
            
    end
end
%%
% target enhnaced cells for each target, catch experiments only
val = 5;
inv = 6;
nc = 0;
resp(hits).enhance = cell(1,length(c_allDirs));
resp(fas).enhance = cell(1,length(c_allDirs));
resp(val).enhance = cell(1,length(c_allDirs));
resp(inv).enhance = cell(1,length(c_allDirs));
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        if mouse(imouse).expt(iexp).info.isCatch
        ind = find(ismember(c_allDirs,chop(mouse(imouse).expt(iexp).catchTargets,2)));
        mod_ind = find(ismember(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).dirs,chop(mouse(imouse).expt(iexp).catchTargets,2)));
        c = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs(mod_ind);
        for idir = 1:length(ind)
            resp(hits).enhance{ind(idir)} = cat(2, resp(hits).enhance{ind(idir)}, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir}(:,c{idir},:),3));
            resp(val).enhance{ind(idir)} = cat(2, resp(val).enhance{ind(idir)}, mean( cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir}(:,c{idir},:),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir}(:,c{idir},:)),3));
            resp(fas).enhance{ind(idir)} = cat(2, resp(fas).enhance{ind(idir)}, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir}(:,c{idir},:),3));
            resp(inv).enhance{ind(idir)} = cat(2, resp(inv).enhance{ind(idir)}, mean( cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir}(:,c{idir},:),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir}(:,c{idir},:)),3));
            
        end
        end
    end
end

h1_tc = cellfun(@(x) bsxfun(@minus,nanmean(x,2),nanmean(nanmean(x(pre_win,:),2),1)),resp(val).enhance,'unif',false);
h1_ste = cellfun(@(x) std(x,[],2)/sqrt(length(x)),resp(hits).enhance,'unif',false);
h2_tc = cellfun(@(x) bsxfun(@minus,nanmean(x,2),nanmean(nanmean(x(pre_win,:),2),1)),resp(inv).enhance,'unif',false);
h2_ste = cellfun(@(x) std(x,[],2)/sqrt(length(x)),resp(fas).enhance,'unif',false);

h1_colors = brewermap(length(c_allDirs)+1,'Greys');
h2_colors = brewermap(length(c_allDirs)+1,'Blues');

% timecourses for enhanced cells for each target direction
figure;
clear h
for idir = 1:length(c_allDirs)
    subplot(1,length(c_allDirs),idir)
    h1{idir} = shadedErrorBar(ttMs,h1_tc{idir},h1_ste{idir});
    h1{idir}.mainLine.Color = h1_colors(idir+1,:);
    h1{idir}.mainLine.LineWidth = 3;
    h1{idir}.patch.FaceAlpha = .25;
    hold on
    h2{idir} = shadedErrorBar(ttMs,h2_tc{idir},h2_ste{idir});
    h2{idir}.mainLine.Color = h2_colors(idir+1,:);
    h2{idir}.mainLine.LineWidth = 3;
    hold on
    ylim([-0.005 0.03])
    title(['target enhanced - ' num2str(c_allDirs(idir)) ' deg']);
    axis square
end
xlabel('ms')
ylabel('dF/F')


print([fnout '_targetEnhanceCells_eachTarget'],'-dpdf')

% target tuning - hits only
tuning_tc_mean = cellfun(@(x) nanmean(x,2),resp(hits).tuning,'unif',false);
tuning_tc_ste = cellfun(@(x) nanstd(x,[],2)/sqrt(size(x,2)),resp(hits).tuning,'unif',false);

tuning_r_mean = cellfun(@(x) bsxfun(@minus,nanmean(nanmean(x(trans_win,:),2),1),nanmean(nanmean(x(pre_win,:),2),1)),resp(hits).tuning,'unif',false);
tuning_r_ste = cellfun(@(x) nanstd(mean(x(trans_win,:),1),[],2)/sqrt(size(x,2)),resp(hits).tuning,'unif',false);

tuning_titles = {'<25';'>25 & <60';'>60'};
figure;
suptitle('target enhanced cells')
for ibin = 1:3
    subplot(1,3,ibin)
    d = cell2mat(tuning_r_mean(ibin,:));
    ind = find(~isnan(d));
    h = errorbar(allDirs(ind),d(ind),cell2mat(tuning_r_ste(ibin,ind)),'k');
    h.LineStyle = 'none';
    h.Marker = 'o';
    h.Parent.XTick = allDirs;
    xlim([0 100]);
    ylim([-0.01 0.03])
    xlabel('target orientation')
    ylabel('dF/F')
    axis square
    title(tuning_titles{ibin})
end
%%
dirs_colors = flipud(brewermap(length(allDirs)+1,'*Greys'));
% dirs_colors = dirs_colors(length(allDirs)+1:end,:);

figure;
for idir = 1:length(allDirs)
   h = cdfplot(auc_allDirs{idir});
   h.Color = dirs_colors(idir+1,:);
   h.LineWidth = 3;
   hold on
end
xlabel('auc')
ylabel('fraction of cells')
title('target enhancement by direction of ori change')
legend(strread(num2str(allDirs),'%s'))

auc_mean_dirs = cell2mat(cellfun(@(x) mean(x,1),auc_allDirs,'unif',false));
auc_ste_dirs = cell2mat(cellfun(@(x) std(x)/sqrt(length(x)),auc_allDirs,'unif',false));

% mean auc for each target direction
figure;
errorbar(allDirs,auc_mean_dirs,auc_ste_dirs,'k')
xlabel('orientation')
ylabel('mean auc')
title('mean auc for each target ori')

prctl = 10;
[auc_dirs_sort auc_dirs_sortInd] = cellfun(@(x) sort(x),auc_allDirs,'unif',false);

topAUC_dirs_ind = cellfun(@(x) x(length(x)-floor(length(x)/prctl):end),auc_dirs_sortInd,'unif',false);

tc_valHits_topAUC_dirs = cell2mat(cellfun(@(x,y) mean(x(:,y),2)-mean(mean(x(pre_win,y),2),1),resp(hits).all,topAUC_dirs_ind,'unif',false));

% timecourses for top auc cells for each target direction
figure;
suptitle(['tc of top ' num2str(prctl) 'prct AUC cells for each target direction']);
clear h
for idir = 1:length(allDirs)
    subplot(2,5,idir)
    h = plot(ttMs,tc_valHits_topAUC_dirs(:,idir));
    h.Color = 'k';%dirs_colors(idir+1,:);
    h.LineWidth = 3;
    hold on
    xlim([-1000 2000])
    ylim([-0.05 0.05])
    xlabel('ms')
    ylabel('dF/F')
    title(num2str(allDirs(idir)))
    axis square
end

%time-courses of example top cells
figure;
suptitle('example top cells')
for idir = 1:length(allDirs)
   subplot(2,5,idir)
   temp = bsxfun(@minus,resp(hits).all{idir},mean(resp(hits).all{idir}(pre_win,:),1));
   ind = randperm(length(topAUC_dirs_ind{idir}),1);
   h = plot(ttMs,temp(:,topAUC_dirs_ind{idir}(ind)),'k');
   h.LineWidth = 3;
   hold on
   vline(0,'k:')
   vline([trans_win(1) trans_win(end)]*(cycTimeMs/cycTime),'r--')
   xlabel('ms')
   ylabel('dF/F')
   xlim([-1000 2000])
   ylim([-0.05 0.1])
   title(num2str(allDirs(idir)))
   axis square
end


% heatmap of all cells sorted by auc for each target direction
ndir = length(allDirs);
nd1 = ceil(sqrt(ndir+1));
if (nd1^2)-nd1 > ndir+1
    nB = nd1-1;
else
    nB = nd1;
end

resp_v_sub = cellfun(@(x) bsxfun(@minus,x,mean(x(pre_win,:),1)),resp(hits).all,'unif',false);

[resp_v_sort resp_v_sortInd] = cellfun(@(x) sort(mean(x(trans_win,:),1)),resp_v_sub,'unif',false);

time_ind = find(ismember(ttMs,[0 1000 2000]));

aucSortFig = figure;
suptitle('resp to each valid target, cells sorted by auc')
colormap(brewermap(1001,'*RdBu'))
respSortFig = figure;
suptitle('resp to each valid target, cells sorted by resp')
colormap(brewermap(1001,'*RdBu'))
for idir = 1:ndir
    figure(aucSortFig);
   subplot(nd1,nB,idir)   
   img = imagesc(resp_v_sub{idir}(1:time_ind(2),auc_dirs_sortInd{idir})');
   img.Parent.XTick = time_ind(1:2);
   img.Parent.XTickLabel = strread(num2str([0 1]),'%s');    
    title(num2str(allDirs(idir)))
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    axis square
    
    figure(respSortFig);
   subplot(nd1,nB,idir)
    img = imagesc(resp_v_sub{idir}(1:time_ind(2),resp_v_sortInd{idir})');
   img.Parent.XTick = time_ind(1:2);
   img.Parent.XTickLabel = strread(num2str([0 1]),'%s');    
    title(num2str(allDirs(idir)))
    xlabel('time (s)')
    ylabel('cells')
    clim_vis = [-.1 0.1];
    caxis(clim_vis)
    colorbar
    axis square
end


%% find target response for 3 bins - targ resp sgnf above, below, and neither from last base stim

nbins = 4;
compInd = 1;

auc_fig = figure;
heatmap_fig = figure;
colormap(brewermap([],'*RdBu'))
for iAuc = 1:nbins
i = 1;

resp_hvsfa_all = [];
resp_favsh_all = [];
resp_mvscr_all = [];
resp_crvsm_all = [];
resp_hvscr_all = [];
resp_crvsh_all = [];
resp_hvsm_all = [];
resp_mvsh_all = [];
resp_favscr_all = [];
resp_crvsfa_all = [];
n_hvsfa = zeros(1,size(expt,2));
n_mvscr = zeros(1,size(expt,2));
n_hvscr = zeros(1,size(expt,2));
n_hvsm = zeros(1,size(expt,2));
n_favscr = zeros(1,size(expt,2));
n_all = zeros(1,size(expt,2));

tc_hvsfa_all = [];
tc_favsh_all = [];
tc_mvscr_all = [];
tc_crvsm_all = [];
tc_hvscr_all = [];
tc_crvsh_all = [];
tc_hvsm_all = [];
tc_mvsh_all = [];
tc_favscr_all = [];
tc_crvsfa_all = [];

resp_val_all = [];
resp_inv_all = [];
resp_all_all = [];
tc_val_all = [];
tc_inv_all = [];
tc_all = [];

%find all catch dirs used
cds = [];
trLs = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)   
        cds = [cds mouse(imouse).expt(iexp).info.cDirs];
        trLs = [trLs mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).tcyc*mouse(imouse).expt(iexp).info.cyc_time_ms];
    end
end
cds = unique(cds);
trLs = unique(trLs);

resp_h_dir = cell(length(cds),1);
resp_fa_dir = cell(length(cds),1);
resp_m_dir = cell(length(cds),1);
resp_cr_dir = cell(length(cds),1);
val_all_dir = cell(length(cds),1);
inv_all_dir = cell(length(cds),1);
n_h_dir = zeros(length(cds),size(expt,2));
n_m_dir = zeros(length(cds),size(expt,2));

oriTuningResp_avg = [];
oriTuningResp_tc = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            val_temp = [];
            inv_temp = [];
            all_temp = [];
            % divide up cells by percentile
            mi_val = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value;
            cell_ind = find(mi_val >= percentile_ind(iAuc) & mi_val < percentile_ind(iAuc+1));
%             % divide up cells by target enhancement
%                 c = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).rst_dirs')));
%             if iAuc == 1 %cells that had bigger response to one of the targets
%                 A = cellfun(@(x) find(x > 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value_dirs,'unif',false);
%                 cell_ind = intersect(c,unique(cat(1,cell2mat(A'))));
%             elseif iAuc == 2 %cells that had smaller response to one of the targets
%                 B = cellfun(@(x) find(x < 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value_dirs,'unif',false);
%                 cell_ind = intersect(c,unique(cat(1,cell2mat(B'))));
%             elseif iAuc == 3 %cells that had equivalent response to one of the targets               
%                 nc = length(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value);
%                 cell_ind = setdiff(1:nc,c);
%             end
        % find number of each trial outcome type
            cdirs = mouse(imouse).expt(iexp).info.cDirs;
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false)),3,length(cdirs));
            nDirs_h = sz(3,:);
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false)),3,length(cdirs));
            nDirs_m = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_m = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp,'unif',false)),3,length(cdirs));
            nDirs_fa = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_fa = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp,'unif',false)),3,length(cdirs));
            nDirs_cr = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_cr = nT(ind3);
            end
            
            if sum(nDirs_fa > 1) > 0 & sum(nDirs_h > 1) > 0  
                mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
                dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
                load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
                oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
                oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
            elseif sum(nDirs_m > 1) > 0 & sum(nDirs_cr > 1) > 0 
                mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
                dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
                load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
                oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
                oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
            end
                
        % fa vs hits matched
            ind = find(lt(nDirs_fa,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvsfa = nDirs_fa;
            elseif isempty(ind)
                nDirs_hvsfa = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_fa(ind) nDirs_h(ind2)]
                nDirs_hvsfa = nT(ind3);
            end
            resp_h_dir_temp = cell(length(cdirs),1);
            resp_fa_dir_temp = cell(length(cdirs),1);
            if any(nDirs_hvsfa > 1)
%             resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa),100);
%             resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa),100);           
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            for iboot = 1:100
                resp_rTh = [];
                resp_rTfa = [];
                rH_dir_temp = cell(length(cdirs),1);
                rFA_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rH_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_h_dir_temp{idir} = cat(4,resp_h_dir_temp{idir},rH_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rFA_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_fa_dir_temp{idir} = cat(4,resp_fa_dir_temp{idir},rFA_dir_temp{idir});
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp(fas).all(:,:,:,iboot) = resp_rTfa;
%             if iboot == 1
%                tc_hvsfa = resp_rTh;
%                tc_favsh = resp_rTfa;
%             end
            end
            resp_h = nanmean(resp_h,4);
            resp(fas).all = nanmean(resp(fas).all,4);
            
%             tc_hvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
%             tc_favsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
            
            
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                ind = find(cds == cdirs(idir));
                resp_h_dir{ind,:} = cat(2,resp_h_dir{ind,:},squeeze(nanmean(nanmean(resp_h_dir_temp{idir},4),3)));
                resp_fa_dir{ind,:} = cat(2,resp_fa_dir{ind,:},squeeze(nanmean(nanmean(resp_fa_dir_temp{idir},4),3)));
                end
            end
            else
            resp_h = [];
            resp(fas).all = [];
%             tc_hvsfa = [];
%             tc_favsh = [];
            end
            
            
            tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(resp_h,3));%tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(tc_hvsfa,3)); 
            tc_favsh_all = cat(2,tc_favsh_all,mean(resp(fas).all,3));%tc_favsh_all = cat(2,tc_favsh_all,mean(tc_favsh,3));  
            
            if any(nDirs_hvsfa > 1)
            resp_hvsfa = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
            resp_favsh = squeeze(mean(mean(resp(fas).all(trans_win,:,:),3),1))-squeeze(mean(mean(resp(fas).all(pre_win,:,:),3),1));
            val_temp = resp_h;
            inv_temp = resp(fas).all;
            all_temp = cat(3,resp_h,resp(fas).all);
            else
                resp_hvsfa = [];
                resp_favsh = [];
                val_temp = [];
                inv_temp = [];
                all_temp = [];
            end

%             resp_val_all = cat(2,resp_val_all,resp_hvsfa);
%             resp_inv_all = cat(2,resp_inv_all,resp_favsh);
%             tc_val_all = cat(2,tc_val_all,mean(resp_h,3));%tc_val_all = cat(2,tc_val_all,mean(tc_hvsfa,3));
%             tc_inv_all = cat(2,tc_inv_all,mean(resp(fas).all,3));%tc_inv_all = cat(2,tc_inv_all,mean(tc_favsh,3));
%             
%             resp_hvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_favsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
            resp_hvsfa_all = cat(2,resp_hvsfa_all,resp_hvsfa);
            resp_favsh_all = cat(2,resp_favsh_all,resp_favsh);
            n_hvsfa(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1));
            for idir = find(nDirs_hvsfa > 1)
                ind = find(cds == cdirs(idir));
                n_h_dir(ind,i) = nDirs_hvsfa(idir);
            end

            ind = find(lt(nDirs_cr,nDirs_m));
            if length(ind) == length(cdirs)
                nDirs_mvscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_mvscr = nDirs_m;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_m(ind2)];
                nDirs_mvscr = nT(ind3);
            end           
            resp_m_dir_temp = cell(length(cdirs),1);
            resp_cr_dir_temp = cell(length(cdirs),1);
            if any(nDirs_mvscr > 1)
            resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
%             resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr),100);
%             resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr),100);
            for iboot = 1:100
                resp_rTm = [];
                resp_rTcr = [];
                rM_dir_temp = cell(length(cdirs),1);
                rCR_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rM_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_m_dir_temp{idir} = cat(4,resp_m_dir_temp{idir},rM_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rCR_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_cr_dir_temp{idir} = cat(4,resp_cr_dir_temp{idir},rCR_dir_temp{idir});
                end
            end
            end
            resp_m(:,:,:,iboot) = resp_rTm;
            resp_cr(:,:,:,iboot) = resp_rTcr;            
%             if iboot == 1
%                tc_mvscr = resp_rTm;
%                tc_crvsm = resp_rTcr;
%             end
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                ind = find(cds == cdirs(idir));
                resp_m_dir{ind,:} = cat(2,resp_m_dir{ind,:},squeeze(nanmean(nanmean(resp_m_dir_temp{idir},4),3)));
                resp_cr_dir{ind,:} = cat(2,resp_cr_dir{ind,:},squeeze(nanmean(nanmean(resp_cr_dir_temp{idir},4),3)));
                end
            end
            
            resp_m = nanmean(resp_m,4);
            resp_cr = nanmean(resp_cr,4);
            else
                resp_m = [];
                resp_cr = [];
%                 tc_mvscr = [];
%                 tc_crvsm = [];
            end
            
            tc_mvscr_all = cat(2,tc_mvscr_all,mean(resp_m,3)); 
            tc_crvsm_all = cat(2,tc_crvsm_all,mean(resp_cr,3)); 
            
            if any(nDirs_mvscr > 1)
                    resp_mvscr = squeeze(mean(mean(resp_m(trans_win,:,:),3),1))-squeeze(mean(mean(resp_m(pre_win,:,:),3),1));
                resp_crvsm = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
                val_temp = cat(3,val_temp,resp_m);
                inv_temp = cat(3,inv_temp,resp_cr);
                all_temp = cat(3,all_temp,resp_m,resp_cr);
            else
                resp_mvscr = [];
                resp_crvsm = [];
                val_temp = val_temp;
                inv_temp = inv_temp;
                all_temp = all_temp;
            end
            
            if ~isempty(val_temp) | ~isempty(inv_temp)
            resp_val_all = cat(2,resp_val_all,squeeze(mean(mean(val_temp(trans_win,:,:),3),1))-squeeze(mean(mean(val_temp(pre_win,:,:),3),1)));%cat(2,resp_val_all,resp_mvscr);
            resp_inv_all = cat(2,resp_inv_all,squeeze(mean(mean(inv_temp(trans_win,:,:),3),1))-squeeze(mean(mean(inv_temp(pre_win,:,:),3),1)));%cat(2,resp_inv_all,resp_crvsm);
            resp_all_all = cat(2,resp_all_all,squeeze(mean(mean(all_temp(trans_win,:,:),3),1))-squeeze(mean(mean(all_temp(pre_win,:,:),3),1)));
            tc_val_all = cat(2,tc_val_all,mean(val_temp,3));%mean(resp_m,3));
            tc_inv_all = cat(2,tc_inv_all,mean(inv_temp,3));%mean(resp_cr,3));
            tc_all = cat(2,tc_all,mean(all_temp,3));
            end
            n_all(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1))+sum(nDirs_mvscr(nDirs_mvscr > 1));
            
%             resp_mvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_mvscr_all = cat(2,resp_mvscr_all,resp_mvscr);
            resp_crvsm_all = cat(2,resp_crvsm_all,resp_crvsm);
            n_mvscr(i) = sum(nDirs_mvscr(nDirs_mvscr > 1));
            for idir = find(nDirs_mvscr > 1)
                ind = find(cds == cdirs(idir));
                n_m_dir(ind,i) = nDirs_mvscr(idir);
            end

        tH = cellfun(@(resp_h_dir_temp) nanmean(resp_h_dir_temp,4),resp_h_dir_temp,'unif',false);
        tM = cellfun(@(resp_m_dir_temp) nanmean(resp_m_dir_temp,4),resp_m_dir_temp,'unif',false);
        tFA = cellfun(@(resp_fa_dir_temp) nanmean(resp_fa_dir_temp,4),resp_fa_dir_temp,'unif',false);
        tCR = cellfun(@(resp_cr_dir_temp) nanmean(resp_cr_dir_temp,4),resp_cr_dir_temp,'unif',false);
        
        for idir = 1:length(cdirs)
            tVal = mean(cat(3,tH{idir},tM{idir}),3);
            tInv = mean(cat(3,tFA{idir},tCR{idir}),3);
            ind = find(cds == cdirs(idir));
            val_all_dir{ind} = cat(2,val_all_dir{ind},tVal);
            inv_all_dir{ind} = cat(2,inv_all_dir{ind},tInv);
        end

        % cr vs hits matched
%             hits = 2;
%             crs = 8;

            ind = find(lt(nDirs_cr,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_hvscr = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_h(ind2)];
                nDirs_hvscr = nT(ind3);
            end
            if any(nDirs_hvscr > 1)
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)),100);
%             tc_hvscr = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
%             tc_crvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
            for iboot = 1:100
            resp_rTh = [];
            resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_hvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_cr(:,:,:,iboot) = resp_rTcr;
%             if iboot == 1
%                tc_hvscr = resp_rTh;
%                tc_crvsh = resp_rTcr;
%             end
            end
            resp_h = mean(resp_h,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_h = [];
            resp_cr = [];
%             tc_hvscr = [];
%             tc_crvsh = [];
            end   
            
            tc_hvscr_all = cat(2,tc_hvscr_all,mean(resp_h,3)); 
            tc_crvsh_all = cat(2,tc_crvsh_all,mean(resp_cr,3));  
            
            
            if any(nDirs_hvscr > 1)
                resp_hvscr = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
                resp_crvsh = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_hvscr = [];
                resp_crvsh = [];
            end

%             resp_hvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_hvscr_all = cat(2,resp_hvscr_all,resp_hvscr);
            resp_crvsh_all = cat(2,resp_crvsh_all,resp_crvsh);
            n_hvscr(i) = sum(nDirs_hvscr(nDirs_hvscr > 1));

            ind = find(lt(nDirs_m,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvsm = nDirs_m;
            elseif isempty(ind)
                nDirs_hvsm = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_m(ind) nDirs_h(ind2)];
                nDirs_hvsm = nT(ind3);
            end
            if any(nDirs_hvsm > 1)
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)),100);
            resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)),100);
%             tc_hvsm = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
%             tc_mvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
            for iboot = 1:100
                resp_rTh = [];
                resp_rTm = [];
            for idir = 1:length(cdirs)
                if nDirs_hvsm(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_m(:,:,:,iboot) = resp_rTm;
%             if iboot == 1
%                tc_hvsm = resp_rTh;
%                tc_mvsh = resp_rTm;
%             end
            end
            resp_h = mean(resp_h,4);
            resp_m = mean(resp_m,4);
            else
            resp_h = [];
            resp_m = [];
%             tc_hvsm = [];
%             tc_mvsh = [];
            end
            
            
            tc_hvsm_all = cat(2,tc_hvsm_all,mean(resp_h,3)); 
            tc_mvsh_all = cat(2,tc_mvsh_all,mean(resp_m,3)); 
            
            
            if any(nDirs_hvsm > 1)
                resp_hvsm = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
                resp_mvsh = squeeze(mean(mean(resp_m(trans_win,:,:),3),1))-squeeze(mean(mean(resp_m(pre_win,:,:),3),1));
            else
                resp_hvsm = [];
                resp_mvsh = [];
            end

%             resp_hvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_mvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
            resp_hvsm_all = cat(2,resp_hvsm_all,resp_hvsm);
            resp_mvsh_all = cat(2,resp_mvsh_all,resp_mvsh);
            n_hvsm(i) = sum(nDirs_hvsm(nDirs_hvsm > 1));

            ind = find(lt(nDirs_cr,nDirs_fa));
            if length(ind) == length(cdirs)
                nDirs_favscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_favscr = nDirs_fa;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_fa(ind2)];
                nDirs_favscr = nT(ind3);
            end
            if any(nDirs_favscr > 1)
            resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);           
%             tc_favscr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
%             tc_crvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
            for iboot = 1:100
                resp_rTfa = [];
                resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_favscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                end
            end
            resp(fas).all(:,:,:,iboot) = resp_rTfa;
            resp_cr(:,:,:,iboot) = resp_rTcr;
%             if iboot == 1
%                tc_favscr = resp_rTfa;
%                tc_crvsfa = resp_rTcr;
%             end
            end
            resp(fas).all = mean(resp(fas).all,4);
            resp_cr = mean(resp_cr,4);
            else
            resp(fas).all = [];
            resp_cr = [];
%             tc_favscr = [];
%             tc_crvsfa = [];
            end
            
            tc_favscr_all = cat(2,tc_favscr_all,mean(resp(fas).all,3)); 
            tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(resp_cr,3)); 
            
            
            if any(nDirs_favscr > 1)
                resp_favscr = squeeze(mean(mean(resp(fas).all(trans_win,:,:),3),1))-squeeze(mean(mean(resp(fas).all(pre_win,:,:),3),1));
                resp_crvsfa = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_favscr = [];
                resp_crvsfa = [];
            end

%             resp_favscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_favscr_all = cat(2,resp_favscr_all,resp_favscr);
            resp_crvsfa_all = cat(2,resp_crvsfa_all,resp_crvsfa);
            n_favscr(i) = sum(nDirs_favscr(nDirs_favscr > 1));

          
        end
        i = i+1;
    end
end



%scatter
figure(auc_fig);
[h,p] = ttest(resp_val_all,resp_inv_all,'alpha',0.05/(size(resp_inv_all,2)));
subplot(nbins,3,1+((iAuc-1)*3))
scatter(resp_val_all,resp_inv_all,50,'k.')
hold on
errorbarxy(mean(resp_val_all),mean(resp_inv_all),std(resp_val_all)/length(resp_val_all),std(resp_inv_all)/length(resp_inv_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid all')
ylabel('invalid all')
title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells']})

%cdf
figure(auc_fig);
subplot(nbins,3,2+((iAuc-1)*3))
hVvsInv = cdfplot(resp_val_all);
hVvsInv.Color = 'k';
hold on
hInvvsV = cdfplot(resp_inv_all);
hInvvsV.Color = 'c';
[h p] = kstest2(resp_val_all,resp_inv_all);
xlim([-0.05 0.05])
title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
axis square

%time-courses
figure(auc_fig)
subplot(nbins,3,3+((iAuc-1)*3))
tc1 = mean(tc_val_all,2)-mean(mean(tc_val_all(pre_win,:),2));
err1 = std(tc_val_all,[],2)/sqrt(size(tc_val_all,2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = mean(tc_inv_all,2)-mean(mean(tc_inv_all(pre_win,:),2));
err2 = std(tc_inv_all,[],2)/sqrt(size(tc_inv_all,2));
shadedErrorBar(ttMs,tc2,err2,'c');
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];[num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']});
axis square

figure(heatmap_fig);
subplot(nbins,2,1+((iAuc-1)*2))
temp_hm = bsxfun(@minus,tc_val_all,mean(tc_val_all(pre_win,:),1));
% hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*3000),:)');
% hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*2000)];
% hm.Parent.XTickLabel = [0 1 2]
hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000),:)');
hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*250) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*500)];
hm.Parent.XTickLabel = [0 0.25 0.5]
xlabel('time (s)')
ylabel('cells')
title('valid')
axis square
clim([-0.1 0.1])
colorbar
subplot(nbins,2,2+((iAuc-1)*2))
temp_hm = bsxfun(@minus,tc_inv_all,mean(tc_inv_all(pre_win,:),1));
% hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*3000),:)');
% hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*2000)];
% hm.Parent.XTickLabel = [0 1 2]
hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000),:)');
hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*250) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*500)];
hm.Parent.XTickLabel = [0 0.25 0.5]
xlabel('time (s)')
ylabel('cells')
title('invalid')
axis square
clim([-0.1 0.1])
colorbar

end


print([fnout '_auROCBins_val-inv-all'],'-dpdf')


load(fullfile(rc.caOutputDir,dataGroup, [titleStr '_' mouse_str '_modCells']))

trType = 3;
cellInd_auROC = 14;

set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
% set(0,'DefaultaxesFontSize', 16)

auc_visTarResp = [];            
auc_visTarResp_short = [];
auc_visTarResp_long = [];
rs_p_visTarResp = []; 
rs_visTarResp = []; 
rs_visTarResp_maxDir = [];
auc_visTarRespMaxDir = [];
expColor = [];
i=0;
% rh_base = [];
% rh_tar = [];
targetBaseRatio = [];
tc = [];
base_resp = [];
tar_resp = [];
taskCells = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch 
            if any(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false))) | any(cell2mat(cellfun(@(x) ~isempty(x),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false)))
            %auROC across all trials
            a = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value;           
            as = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(1).value;            
            al = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(2).value;
            u = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst_p;
            ut = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst; 
            ut_dirs = cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst_dirs_logical);
            ut_dirs = logical(ut_dirs(:));
            ut_maxDir_ind = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst_dirs')));
            ut_maxDir = zeros(1,length(a));
            ut_maxDir(ut_maxDir_ind) = 1;
            
            if cellInd_auROC == 15;                
            cell_ind = intersect(mouse(imouse).expt(iexp).cells(8).ind,mouse(imouse).expt(iexp).cells(12).ind);
            else
            cell_ind = mouse(imouse).expt(iexp).cells(cellInd_auROC).ind;
            end
            cell_v = zeros(length(a),1);
            cell_v(cell_ind) = 1;
            taskCells = cat(1,taskCells,cell_v);
            
            % max or min auROC across by each target (biggest effect) *has
            % to be significant rst
%             b = max(cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value_dirs),[],2);
            b_sub =  cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value_dirs)-0.5;
            b_sub_rect = zeros(size(b_sub(:)));
            b_sub_rect(ut_dirs) = b_sub(ut_dirs);
            b_sub_rect = reshape(b_sub_rect,size(b_sub));
            [b_abs b_maxInd] = max(abs(b_sub_rect),[],2);
            b = (b_sub(sub2ind(size(b_sub),1:size(b_sub,1),b_maxInd'))+0.5)';
            
            if size(a,1) > size(a,2)
                a = a';
            end
            auc_visTarResp = cat(2,auc_visTarResp,a);            
            auc_visTarResp_short = cat(1,auc_visTarResp_short,as);
            auc_visTarResp_long = cat(2,auc_visTarResp_long,al);
            rs_visTarResp = cat(2,rs_visTarResp,ut); 
            rs_visTarResp_maxDir = cat(2,rs_visTarResp_maxDir,ut_maxDir); 
            rs_p_visTarResp = cat(2,rs_p_visTarResp,u); 
            auc_visTarRespMaxDir = cat(1,auc_visTarRespMaxDir,b);
            expColor = cat(1,expColor,ones(length(a),1)+i);
            i = i+1;
            
            %resp histogram each cell
%             rh_base = cat(1,rh_base,msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).resp_base);
%             rh_tar = cat(1,rh_tar,msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).resp_tar);
            targetBaseRatio = cat(2,targetBaseRatio,msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).ratio);
            resp = mean(cat(3,cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{:}),cat(3,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{:})),3);
            tc = cat(2,tc,resp);
            base_resp = cat(2,base_resp,mean(resp(trans_win-cycTime,:),1)-mean(resp(pre_win-cycTime,:),1));
            tar_resp = cat(2,tar_resp,mean(resp(trans_win,:),1)-mean(resp(pre_win,:),1));
            end
        end
    end
end

base_resp_rect = base_resp - min(base_resp);
tar_resp_rect = tar_resp - min(tar_resp);
targetBaseRatio_rect = bsxfun(@rdivide,tar_resp_rect,base_resp_rect);

c = length(auc_visTarResp);
nbins = 4;
val_sort = sort(auc_visTarResp);
percentile_ind = [val_sort(1:floor(c/nbins):c) val_sort(end)];
valMax_sort = sort(auc_visTarRespMaxDir);
perc_ind_maxDir = [valMax_sort(1:floor(c/nbins):c)' valMax_sort(end)];

auROCall_fig = figure;
subplot(2,3,1)
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

figure(auROCall_fig);
subplot(2,3,2)
colormap(jet)
scatter(auc_visTarResp,rs_p_visTarResp,50,expColor);
hold on
xlim([0 1])
hline(0.05,'k--')
xlabel('auc')
ylabel('p value - rank sum test')
colorbar
clim([1 8])

%color for signif rank sum test
figure(auROCall_fig);
subplot(2,3,3)
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
title('1 = signif rs test')
axis square

% colorbar for ea exp
figure(auROCall_fig);
subplot(2,3,4)
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

%plot (valid only) timecourses for cells that are signif and those that
%aren't
tc = bsxfun(@minus,tc,mean(tc(pre_win,:),1));
taskCells_ind = find(taskCells);
tc_taskCells = tc(:,taskCells_ind);
ntrCells_ind = setdiff(1:size(tc,2),taskCells_ind);

figure;
suptitle('experiments with invalid trials')
subplot(2,2,1)
h = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5))),2),std(tc(:,intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp),find(auc_visTarRespMaxDir > 0.5)))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5))),2),std(tc(:,intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp),find(auc_visTarRespMaxDir < 0.5)))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,find(rs_visTarResp == 0)),2),std(tc(:,find(rs_visTarResp == 0)),[],2)/sqrt(length(find(rs_visTarResp == 0))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('average auROC across all targets')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5))))], ['+rs < 0.5, ' num2str(length(intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5))))], ['ns, ' num2str(length(find(rs_visTarResp == 0)))]},'Location','NorthWest')

subplot(2,2,2)
h(1) = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5))),2),std(tc(:,intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5)))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5))),2),std(tc(:,intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5)))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,find(rs_visTarResp_maxDir == 0)),2),std(tc(:,find(rs_visTarResp_maxDir == 0)),[],2)/sqrt(length(find(rs_visTarResp_maxDir == 0))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('max auROC across all targets')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5))))], ['+rs < 0.5, ' num2str(length(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5))))], ['ns, ' num2str(length(find(rs_visTarResp_maxDir == 0)))]},'Location','NorthWest')

subplot(2,2,3)
h = shadedErrorBar(ttMs,mean(tc_taskCells(:,intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) > 0.5))),2),std(tc_taskCells(:,intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) > 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) > 0.5)))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc_taskCells(:,intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) < 0.5))),2),std(tc_taskCells(:,intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) < 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) < 0.5)))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc_taskCells(:,find(rs_visTarResp(taskCells_ind) == 0)),2),std(tc_taskCells(:,find(rs_visTarResp(taskCells_ind) == 0)),[],2)/sqrt(length(find(rs_visTarResp(taskCells_ind) == 0))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('task responsive cells - avg')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) > 0.5))))], ['+rs < 0.5, ' num2str(length(intersect(find(rs_visTarResp(taskCells_ind)),find(auc_visTarResp(taskCells_ind) < 0.5))))], ['ns, ' num2str(length(find(rs_visTarResp(taskCells_ind) == 0)))]},'Location','NorthWest')

subplot(2,2,4)
h(1) = shadedErrorBar(ttMs,mean(tc_taskCells(:,intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) > 0.5))),2),std(tc_taskCells(:,intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) > 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) > 0.5)))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc_taskCells(:,intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) < 0.5))),2),std(tc_taskCells(:,intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) < 0.5))),[],2)/sqrt(length(intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) < 0.5)))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc_taskCells(:,find(rs_visTarResp_maxDir(taskCells_ind) == 0)),2),std(tc_taskCells(:,find(rs_visTarResp_maxDir(taskCells_ind) == 0)),[],2)/sqrt(length(find(rs_visTarResp_maxDir(taskCells_ind) == 0))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('task responsive cells - max')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) > 0.5))))], ['+rs < 0.5, ' num2str(length(intersect(find(rs_visTarResp_maxDir(taskCells_ind)),find(auc_visTarRespMaxDir(taskCells_ind) < 0.5))))], ['ns, ' num2str(length(find(rs_visTarResp_maxDir(taskCells_ind) == 0)))]},'Location','NorthWest')

print([fnout '_tcAllCellsByRST'],'-dpdf','-fillpage')

% time-courses of non-task responsive cells
figure
suptitle('non-task responsive cells only')
subplot(1,2,1)
h = shadedErrorBar(ttMs,mean(tc(:,intersect(intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5)),ntrCells_ind)),2),std(tc(:,intersect(intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5)),ntrCells_ind)),[],2)/sqrt(length(intersect(intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5)),ntrCells_ind))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5)),ntrCells_ind)),2),std(tc(:,intersect(intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5)),ntrCells_ind)),[],2)/sqrt(length(intersect(intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5)),ntrCells_ind))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp == 0),ntrCells_ind)),2),std(tc(:,intersect(find(rs_visTarResp == 0),ntrCells_ind)),[],2)/sqrt(length(intersect(find(rs_visTarResp == 0),ntrCells_ind))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('average auROC across all targets')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(intersect(find(rs_visTarResp),find(auc_visTarResp > 0.5)),ntrCells_ind)))], ['+rs < 0.5, ' num2str(length(intersect(intersect(find(rs_visTarResp),find(auc_visTarResp < 0.5)),ntrCells_ind)))], ['ns, ' num2str(length(intersect(find(rs_visTarResp == 0),ntrCells_ind)))]},'Location','NorthWest')

subplot(1,2,2)
h = shadedErrorBar(ttMs,mean(tc(:,intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5)),ntrCells_ind)),2),std(tc(:,intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5)),ntrCells_ind)),[],2)/sqrt(length(intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5)),ntrCells_ind))));
h.mainLine.Color = [0.75 0 0];
h.mainLine.LineWidth = 3;
leg(1) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5)),ntrCells_ind)),2),std(tc(:,intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5)),ntrCells_ind)),[],2)/sqrt(length(intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5)),ntrCells_ind))));
h.mainLine.Color = [0 0 0.75];
h.mainLine.LineWidth = 3;
leg(2) = h.mainLine;
hold on
h = shadedErrorBar(ttMs,mean(tc(:,intersect(find(rs_visTarResp_maxDir == 0),ntrCells_ind)),2),std(tc(:,intersect(find(rs_visTarResp_maxDir == 0),ntrCells_ind)),[],2)/sqrt(length(intersect(find(rs_visTarResp_maxDir == 0),ntrCells_ind))));
h.mainLine.LineWidth = 3;
leg(3) = h.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('max auROC across all targets')
legend(leg,{['+rs > 0.5, ' num2str(length(intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir > 0.5)),ntrCells_ind)))], ['+rs < 0.5, ' num2str(length(intersect(intersect(find(rs_visTarResp_maxDir),find(auc_visTarRespMaxDir < 0.5)),ntrCells_ind)))], ['ns, ' num2str(length(intersect(find(rs_visTarResp_maxDir == 0),ntrCells_ind)))]},'Location','NorthWest')

print([fnout '_tcNonTaskResponsiveCellsByRST'],'-dpdf','-fillpage')

%plot resp ratio scatter for each cell, sorted by auc
respRatioFig = figure;
suptitle({'mean of resp window for target and last base stim'; 'colorbar is auROC';'task responsive cells only'})
subplot(1,2,1)
colormap(brewermap([],'*RdBu'))
s = scatter(base_resp(taskCells_ind),tar_resp(taskCells_ind),30,double(auc_visTarResp(taskCells_ind)),'o','filled');
hold on
ss = scatter(base_resp(intersect(find(rs_visTarResp),taskCells_ind)),tar_resp(intersect(find(rs_visTarResp),taskCells_ind)),30,double(auc_visTarResp(intersect(find(rs_visTarResp),taskCells_ind))),'o','filled');
ss.MarkerEdgeColor = [0 0 0];
hold on
plot([-1:0.1:1],[-1:0.1:1],'k--')
xlim([-0.05 0.12])
ylim([-0.05 0.12])
xlabel('base')
ylabel('target')
colorbar
caxis([0 1])
axis square
title('avg auROC')
subplot(1,2,2)
colormap(brewermap([],'*RdBu'))
s = scatter(base_resp(taskCells_ind),tar_resp(taskCells_ind),30,auc_visTarRespMaxDir(taskCells_ind),'o','filled');
hold on
ss = scatter(base_resp(intersect(find(rs_visTarResp_maxDir),taskCells_ind)),tar_resp(intersect(find(rs_visTarResp_maxDir),taskCells_ind)),30,auc_visTarRespMaxDir(intersect(find(rs_visTarResp_maxDir),taskCells_ind)),'o','filled');
ss.MarkerEdgeColor = [0 0 0];
hold on
plot([-1:0.1:1],[-1:0.1:1],'k--')
l = vline(0,'-');
l.Color = [0.5 0.5 0.5];
l = hline(0,'-');
l.Color = [0.5 0.5 0.5]
xlim([-0.05 0.12])
ylim([-0.05 0.12])
xlabel('base')
ylabel('target')
colorbar
caxis([0 1])
axis square
title('max auROC across all targets')

print([fnout '_ratioScatter_auROCcolor'],'-dpdf','-fillpage')
%
maxAUC_noise = find(auc_visTarRespMaxDir == 1);

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

clear resp
auc_allDirs = cell(1,length(allDirs));
auc_allDirs_catch = cell(1,length(allDirs));
auc_avg_catch = cell(1,length(allDirs));
base_resp_dir = cell(1,length(allDirs));
tar_resp_dir = cell(1,length(allDirs));
rs_dirs = cell(1,length(allDirs));
rs_dirs_avg = cell(1,length(allDirs));
rs_count_val = [];
rs_count_inv = [];
rs_count_perOri = zeros(1,length(allDirs));
rs_count_perOri_rsavg = zeros(1,length(allDirs));
rs_count_perOri_rsmax = zeros(1,length(allDirs));
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
dir_counts1 = zeros(1,length(allDirs));
dir_counts2 = zeros(1,length(allDirs));
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        cell_ind = mouse(imouse).expt(iexp).cells(cellInd_auROC).ind;
        c = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs')));
%         c = intersect(c,cell_ind);
        A = cellfun(@(x) find(x > 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs,'unif',false);
        B = cellfun(@(x) find(x < 0.5),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs,'unif',false);
        ca = intersect(c,unique(cat(1,cell2mat(A'))));
        cb = intersect(c,unique(cat(1,cell2mat(B'))));
        ind = find(ismember(allDirs,msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).dirs));
        rdc = zeros(1,length(ind));
        rdc_rsavg = zeros(1,length(ind));
        rdc_rsmax = zeros(1,length(ind));
        c_all = find(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst);
%         c_all = intersect(c_all,cell_ind);
        for idir = 1:length(ind)
            a = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs{idir}(cell_ind);
            rdc_rsavg(idir) = sum(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs_logical{idir}(c_all));
            rdc_rsmax(idir) = sum(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs_logical{idir}(c));
            rdc(idir) = sum(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs_logical{idir});
            
            auc_allDirs{ind(idir)} = cat(1, auc_allDirs{ind(idir)},a); 
            
            resp(hits).all{ind(idir)} = cat(2, resp(hits).all{ind(idir)}, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp{idir+1}(:,cell_ind,:),3));
%             resp(hits).enhance{ind(idir)} = cat(2, resp(hits).enhance{ind(idir)}, mean(r{idir},3));
        end
            rs_count_perOri(ind) = bsxfun(@plus,rs_count_perOri(ind), rdc);
            rs_count_perOri_rsavg(ind) = bsxfun(@plus,rs_count_perOri_rsavg(ind), rdc_rsavg);
            rs_count_perOri_rsmax(ind) = bsxfun(@plus,rs_count_perOri_rsmax(ind), rdc_rsmax);
        if mouse(imouse).expt(iexp).info.isCatch
            rs_count_val = cat(1,rs_count_val,sum(cat(1,cell2mat(cellfun(@(x) ismember(cell_ind,x),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs,'unif',false))),2));
            if ~isempty(msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).rst_dirs)
                rs_count_inv = cat(1,rs_count_inv,sum(cat(1,cell2mat(cellfun(@(x) ismember(cell_ind,x),msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).rst_dirs,'unif',false))),2));
            end
            for idir = 1:length(ind)
                a = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value_dirs{idir}(cell_ind);
                aa = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).value(cell_ind);
                rsa = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst(cell_ind);
                rs = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs_logical{idir}(cell_ind);
                auc_allDirs_catch{ind(idir)} = cat(1, auc_allDirs_catch{ind(idir)},a);
                auc_avg_catch{ind(idir)} = cat(2, auc_avg_catch{ind(idir)},aa);
                resp_dirs = mean(cat(3,cat(3,mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp{idir}),cat(3,mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).stimResp{idir})),3);
%                 tc = cat(2,tc,resp);
                base_resp_dir{ind(idir)} = cat(2,base_resp_dir{ind(idir)},mean(resp_dirs(trans_win-cycTime,cell_ind),1)-mean(resp_dirs(pre_win-cycTime,cell_ind),1));
                tar_resp_dir{ind(idir)} = cat(2,tar_resp_dir{ind(idir)},mean(resp_dirs(trans_win,cell_ind),1)-mean(resp_dirs(pre_win,cell_ind),1));
                rs_dirs{ind(idir)} = cat(1,rs_dirs{ind(idir)},rs);
                rs_dirs_avg{ind(idir)} = cat(2,rs_dirs_avg{ind(idir)},rsa);
                dir_counts1(ind(idir)) = dir_counts1(ind(idir))+length(a);
                dir_counts2(ind(idir)) = dir_counts2(ind(idir))+length(mean(resp_dirs(trans_win,cell_ind),1));
            end
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
            rsInd = intersect(rsInd,cell_ind);
%             resp(hits).tuning(ibin,ind) = cellfun(@(x,y) cat(2,x,y), resp(hits).tuning(ibin,ind), cellfun(@(z) mean(z(:,rsInd,:),3),mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).stimResp(exInd),'unif',false),'unif',false);
        end
        if mouse(imouse).expt(iexp).info.isCatch
        cind = find(ismember(c_allDirs,msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).dirs));
%         trInd = find(~isempty(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp));
%         r = cellfun(@(x,y) x(:,y,:),mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp(trInd),msModCells(imouse).expt(iexp).av(iav).outcome(fas).mi(1).comp(1).auc(3).rst_dirs,'unif',false);
        for idir = 1:length(cind)
            resp(fas).all{cind(idir)} = cat(2, resp(fas).all{cind(idir)}, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir}(:,cell_ind,:),3));
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

%histogram of number of targets each cell has a significant auROC
bins = 1:length(allDirs);
bin_vals = histc(rs_count_val,bins);
bin_invs = histc(rs_count_inv,bins);

figure(auROCall_fig)
subplot(2,3,5)
bar_rsCounts = bar(bins,cat(2,bin_vals, bin_invs));
bar_rsCounts(1).FaceColor = 'k';
bar_rsCounts(2).FaceColor = 'c';
xlabel('n of signif targets')
ylabel('n cells')
title('n oris with signif auROC for each cell')
legend({'valid','invalid'})
axis square

print([fnout '_auROC_allcells'],'-dpdf','-fillpage')

figure
bar_rsCounts = bar(1:length(allDirs),cat(2,rs_count_perOri',rs_count_perOri_rsavg',rs_count_perOri_rsmax'),1);
hold on
cellfun(@(x) set(x,'XTickLabel',allDirs),{bar_rsCounts(:).Parent},'unif',false)
cellfun(@(x) set(x,'EdgeColor',[1 1 1]),{bar_rsCounts(:)},'unif',false)
bar_rsCounts(1).FaceColor = 'k';
bar_rsCounts(2).FaceColor = [1 0.5 1];
bar_rsCounts(3).FaceColor = [0.75 0 0.75];
xlabel('ori change')
ylabel('n cells')
title({'n cells with signif auROC for each ori';'blk=all;ltpnk=rst signif(avg);dkpnk=rst signif(max)'})
axis square
print([fnout '_rstCellsPerOri'],'-dpdf','-fillpage')

% scatters of responses for each target direction (catch experiments only)
ndir = length(allDirs);
nd1 = ceil(sqrt(ndir+1));
if (nd1^2)-nd1 > ndir+1
    nB = nd1-1;
else
    nB = nd1;
end

targetBaseRatioDirFig_max = figure;
targetBaseRatioDirFig_avg = figure;
suptitle({'mean of resp window for target and last base stim'; 'colorbar is auROC';'task responsive cells onlyl'})
colormap(brewermap([],'*RdBu'));
for i = 1:ndir
    figure(targetBaseRatioDirFig_max);
    subplot(nd1,nB,i)
    s = scatter(base_resp_dir{i},tar_resp_dir{i},30,auc_allDirs_catch{i},'o','filled');
    hold on
    ss = scatter(base_resp_dir{i}(logical(rs_dirs{i})),tar_resp_dir{i}(logical(rs_dirs{i})),30,auc_allDirs_catch{i}(logical(rs_dirs{i})),'o','filled');
    ss.MarkerEdgeColor = [0 0 0];
    hold on
    plot([-1:0.1:1],[-1:0.1:1],'k--')
    xlim([-0.05 0.12])
    ylim([-0.05 0.12])
    xlabel('base')
    ylabel('target')
    colorbar
    caxis([0 1])
    axis square
    title(['max auROC - ' num2str(allDirs(i)) ' targets'])
    
    figure(targetBaseRatioDirFig_avg);
    subplot(nd1,nB,i)
    s = scatter(base_resp_dir{i},tar_resp_dir{i},30,auc_avg_catch{i},'o','filled');
    hold on
    ss = scatter(base_resp_dir{i}(logical(rs_dirs_avg{i})),tar_resp_dir{i}(logical(rs_dirs_avg{i})),30,auc_avg_catch{i}(logical(rs_dirs_avg{i})),'o','filled');
    ss.MarkerEdgeColor = [0 0 0];
    hold on
    plot([-1:0.1:1],[-1:0.1:1],'k--')
    xlim([-0.05 0.12])
    ylim([-0.05 0.12])
    xlabel('base')
    ylabel('target')
    colorbar
    caxis([0 1])
    axis square
    title(['avg auROC - ' num2str(allDirs(i)) ' targets'])
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
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        ind = find(ismember(c_allDirs,chop(mouse(imouse).expt(iexp).catchTargets,2)));
        mod_ind = find(ismember(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).dirs,chop(mouse(imouse).expt(iexp).catchTargets,2)));
        c = cellfun(@(x) intersect(x,cell_ind),msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(3).rst_dirs(mod_ind),'unif',false);
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
clear h1 h2
for idir = 1:length(c_allDirs)
    subplot(1,length(c_allDirs),idir)
    if ~isnan(mean(h1_tc{idir}))
    h1{idir} = shadedErrorBar(ttMs,h1_tc{idir},h1_ste{idir});
    h1{idir}.mainLine.Color = h1_colors(idir+1,:);
    h1{idir}.mainLine.LineWidth = 3;
    h1{idir}.patch.FaceAlpha = .25;
    hold on
    h2{idir} = shadedErrorBar(ttMs,h2_tc{idir},h2_ste{idir});
    h2{idir}.mainLine.Color = h2_colors(idir+1,:);
    h2{idir}.mainLine.LineWidth = 3;
    end
    hold on
    ylim([-0.005 0.06])
    title(['target enhanced - ' num2str(c_allDirs(idir)) ' deg']);
    axis square
end
xlabel('ms')
ylabel('dF/F')


print([fnout '_targetEnhanceCells_eachTarget'],'-dpdf')

% % target tuning - hits only
% % tuning_tc_mean = cellfun(@(x) nanmean(x,2),resp(hits).tuning,'unif',false);
% % tuning_tc_ste = cellfun(@(x) nanstd(x,[],2)/sqrt(size(x,2)),resp(hits).tuning,'unif',false);
% % 
% % tuning_r_mean = cellfun(@(x) bsxfun(@minus,nanmean(nanmean(x(trans_win,:),2),1),nanmean(nanmean(x(pre_win,:),2),1)),resp(hits).tuning,'unif',false);
% % tuning_r_ste = cellfun(@(x) nanstd(mean(x(trans_win,:),1),[],2)/sqrt(size(x,2)),resp(hits).tuning,'unif',false);
% % 
% % tuning_titles = {'<25';'>25 & <60';'>60'};
% % figure;
% % suptitle('target enhanced cells')
% % for ibin = 1:3
% %     subplot(1,3,ibin)
% %     d = cell2mat(tuning_r_mean(ibin,:));
% %     ind = find(~isnan(d));
% %     h = errorbar(allDirs(ind),d(ind),cell2mat(tuning_r_ste(ibin,ind)),'k');
% %     h.LineStyle = 'none';
% %     h.Marker = 'o';
% %     h.Parent.XTick = allDirs;
% %     xlim([0 100]);
% %     ylim([-0.01 0.03])
% %     xlabel('target orientation')
% %     ylabel('dF/F')
% %     axis square
% %     title(tuning_titles{ibin})
% % end
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
   vline([-cycTimeMs 0],'k:')
   vline(([trans_win(1) trans_win(end)]-pre_event_frames)*(cycTimeMs/cycTime),'r--')
   vline((([trans_win(1) trans_win(end)]-pre_event_frames)*(cycTimeMs/cycTime))-cycTimeMs,'b--');
   xlabel('ms')
   ylabel('dF/F')
   xlim([-1000 1000])
   ylim([-0.05 0.1])
   title(num2str(allDirs(idir)))
   axis square
end

print([fnout '_exCellsHighAuROCeaDir'],'-dpdf','-fillpage')



% heatmap of all cells sorted by auc for each target direction

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

nbins = 3;
compInd = 1;
aucColors = brewermap(nbins,'*RdBu');
tc_val_lastBaseStim = [];
tc_inv_lastBaseStim = [];
dscCellsCount = 0;
start = 1;

auc_fig = figure;
colormap(brewermap([],'*RdBu'));
auc_fig_press_early = figure;
colormap(brewermap([],'*RdBu'));
auc_fig_press_late = figure;
colormap(brewermap([],'*RdBu'));
heatmap_fig = figure;
colormap(brewermap([],'*RdBu'));
valOverlay = figure;
colormap(brewermap([],'*RdBu'));
clear legData;
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

%pre-allocate press align tcs

tc_vis_late = [];
tc_aud_late = [];
tc_vis_early = [];
tc_aud_early = [];

oriTuningResp_avg = [];
oriTuningResp_tc = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            val_temp = [];
            inv_temp = [];
            all_temp = [];
%             % divide up cells by percentile
%             mi_val = msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value;
%             cell_ind = intersect(mouse(imouse).expt(iexp).cells(cellsInd).ind,find(mi_val >= percentile_ind(iAuc) & mi_val < percentile_ind(iAuc+1)));
%             % divide up cells by percentile from max auc each target
%             mi_val = max(cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value_dirs),[],2);
%             maxAUCexp_noise = find(mi_val == 1);
%             cell_ind = intersect(mouse(imouse).expt(iexp).cells(cellsInd).ind,find(mi_val >= perc_ind_maxDir(iAuc) & mi_val < perc_ind_maxDir(iAuc+1)));
% %             cell_ind = find(mi_val >= perc_ind_maxDir(iAuc) & mi_val < perc_ind_maxDir(iAuc+1));
%             cell_ind = setdiff(cell_ind,maxAUCexp_noise);
            % divide up cells by target enhancement
%                 c = unique(cat(1,cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).rst_dirs'))); % max 
                c = find(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).rst); % avg 
%                 ar = cat(1,mouse(imouse).expt(iexp).cells(8).ind,mouse(imouse).expt(iexp).cells(12).ind);% anti responsive
%                 ut_dirs = cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).rst_dirs_logical);
%                 ut_dirs = logical(ut_dirs(:));
%                 auROCmax_sub =  cell2mat(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(1).auc(trType).value_dirs)-0.5;
%                 auROCmax_sub_rect = zeros(size(auROCmax_sub(:)));
%                 auROCmax_sub_rect(ut_dirs) = auROCmax_sub(ut_dirs);
%                 auROCmax_sub_rect = reshape(auROCmax_sub_rect,size(auROCmax_sub));
%                 [auROCmax_abs auROCmax_maxInd] = max(abs(auROCmax_sub_rect),[],2);
%                 auROCmax = (auROCmax_sub(sub2ind(size(auROCmax_sub),1:size(auROCmax_sub,1),auROCmax_maxInd'))+0.5)';
            if iAuc == 3 %cells that had bigger response to one of the targets             
                
                A = find(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value > 0.5);
%                 A = find(auROCmax > 0.5);
                cell_ind = intersect(c,A);
%                 cell_ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(cellInd_auROC).ind);
            elseif iAuc == 1 %cells that had smaller response to one of the targets
                B = find(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value < 0.5);
%                 B = find(auROCmax < 0.5);
                cell_ind = intersect(c,B);
%                 cell_ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(cellInd_auROC).ind);
            elseif iAuc == 2 %cells that had equivalent response to one of the targets               
                nc = length(msModCells(imouse).expt(iexp).av(iav).outcome(hits).mi(1).comp(compInd).auc(3).value);
                cell_ind = setdiff(1:nc,c);
%                 cell_ind = intersect(cell_ind,mouse(imouse).expt(iexp).cells(cellInd_auROC).ind);
            end
        % find number of each trial outcome type
            cdirs = mouse(imouse).expt(iexp).info.cDirs;
%             sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false)),3,length(cdirs));
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false)),3,length(cdirs));
            nDirs_h = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_h = nT(ind3);
            end
%             nDirs_h = sz(3,:);
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
            if iAuc == 3
            dscCellsCount = dscCellsCount+length(cell_ind);
            disp([mouse(imouse).expt(iexp).date '-h' num2str(length(mouse(imouse).expt(iexp).cells(cellInd_auROC).ind))])
            end
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
            if iAuc == 3 & any(nDirs_hvsfa > 1)
                dscCellsCount = dscCellsCount+length(cell_ind);
                disp([mouse(imouse).expt(iexp).date '-m' num2str(length(mouse(imouse).expt(iexp).cells(cellInd_auROC).ind))])
            end
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
            
%********* % press align tc
            cycFrames = mouse(imouse).expt(iexp).info.cyc_time;
            fr_temp = (cycFrames/mouse(imouse).expt(iexp).info.cyc_time_ms);
            minLongMs = 2800;
            minShortMs = 1400;
            minLongFrames = floor(minLongMs*fr_temp);
            minShortFrames = floor(minShortMs*fr_temp);
            cycsLong = ceil(minLongFrames/cycFrames);
            cycsShort = ceil(minShortFrames/cycFrames);
            
            if size(mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp,2) >= cycsLong
            v_long_temp = mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp{cycsLong};
            a_long_temp = mouse(imouse).expt(iexp).align(1).av(2).outcome(hits).cmlvCycResp{cycsLong};
            
            v_sht_temp = mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp{cycsShort};
            a_sht_temp = mouse(imouse).expt(iexp).align(1).av(2).outcome(hits).cmlvCycResp{cycsShort};
            
            tc_vis_late = cat(2,tc_vis_late,mean(v_long_temp(1:minLongFrames+pre_event_frames,cell_ind,:),3));
            tc_aud_late = cat(2,tc_aud_late,mean(a_long_temp(1:minLongFrames+pre_event_frames,cell_ind,:),3));
            tc_vis_early = cat(2,tc_vis_early,mean(v_sht_temp(1:minShortFrames+pre_event_frames,cell_ind,:),3));
            tc_aud_early = cat(2,tc_aud_early,mean(a_sht_temp(1:minShortFrames+pre_event_frames,cell_ind,:),3));
            
            end
       

          
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
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells']})
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells']})

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
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
axis square

%time-courses
    %create tc across all trials, aligned to previous base stim
    nCells = size(tc_val_all,2);
    tc_val_lastBaseAlign(:,start:(start+nCells)-1) = tc_val_all;
    tc_inv_lastBaseAlign(:,start:(start+nCells)-1) = tc_inv_all;
    start = start+nCells;

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
xlim([-15 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
h = vline(([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs))-cycTimeMs,'--r');
h(1).Color = [0.5 0 0];
h(2).Color = [0.5 0 0];
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
hold on
vline(([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs))-cycTimeMs, '--k')
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];[num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']});
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];[num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']});
axis square


figure(valOverlay)
subplot(1,2,1)
hold on
l = shadedErrorBar(ttMs,tc1,err1);
l.mainLine.Color = aucColors(iAuc,:);
l.mainLine.LineWidth = 3;
legData(iAuc) = l.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('valid trials')
legendInfo{iAuc} = ['Bin ' num2str(iAuc)];
subplot(1,2,2)
hold on
l = shadedErrorBar(ttMs,tc2,err2);
l.mainLine.Color = aucColors(iAuc,:);
l.mainLine.LineWidth = 3;
legData(iAuc) = l.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('invalid trials')
legendInfo{iAuc} = ['Bin ' num2str(iAuc)];

figure(heatmap_fig);
subplot(nbins,2,1+((iAuc-1)*2))
temp_hm = bsxfun(@minus,tc_val_all,mean(tc_val_all(pre_win,:),1));
% hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*3000),:)');
% hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*2000)];
% hm.Parent.XTickLabel = [0 1 2]
hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000),:)');
hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*250) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*500)];
hm.Parent.XTickLabel = [0 0.25 0.5];
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
hm.Parent.XTickLabel = [0 0.25 0.5];
xlabel('time (s)')
ylabel('cells')
title('invalid')
axis square
clim([-0.1 0.1])
colorbar

%*********%press align figures
early_win = (11:33)+trans_win(1);
late_win = ((66:88)+pre_event_frames)-(trans_win(1)-pre_event_frames);
resp_v_late = mean(mean(tc_vis_late(late_win,:),1),3);
resp_a_late = mean(mean(tc_aud_late(late_win,:),1),3);
resp_v_early = mean(mean(tc_vis_early(early_win,:),1),3);
resp_a_early = mean(mean(tc_aud_early(early_win,:),1),3);


%scatter
%early analysis window
figure(auc_fig_press_early)
[h,p] = ttest(resp_v_early,resp_a_early,'alpha',0.05/(size(resp_a_early,2)));
subplot(nbins,3,1+((iAuc-1)*3))
scatter(resp_v_early,resp_a_early,50,'k.')
hold on
errorbarxy(mean(resp_v_early),mean(resp_a_early),std(resp_v_early)/length(resp_v_early),std(resp_a_early)/length(resp_a_early),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('visual')
ylabel('auditory')
title({['ROC bin ' num2str(iAuc)];['early-press align; p=' num2str(p) ', ' num2str(size(resp_v_early,2)) '  cells']})

%late analyssi window
figure(auc_fig_press_late)
[h,p] = ttest(resp_v_late,resp_a_late,'alpha',0.05/(size(resp_a_late,2)));
subplot(nbins,3,1+((iAuc-1)*3))
scatter(resp_v_late,resp_a_late,50,'k.')
hold on
errorbarxy(mean(resp_v_late),mean(resp_a_late),std(resp_v_late)/length(resp_v_late),std(resp_a_late)/length(resp_a_late),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('visual')
ylabel('auditory')
title({['ROC bin ' num2str(iAuc)];['late-press align; p=' num2str(p) ', ' num2str(size(resp_v_late,2)) '  cells']})

%cdf
%early
figure(auc_fig_press_early);
subplot(nbins,3,2+((iAuc-1)*3))
hVvsInv = cdfplot(resp_v_early);
hVvsInv.Color = 'g';
hold on
hInvvsV = cdfplot(resp_a_early);
hInvvsV.Color = 'k';
[h p] = kstest2(resp_v_early,resp_a_early);
xlim([-0.05 0.05])
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
title({['ROC bin ' num2str(iAuc)];['early; p= ' num2str(p)]})
axis square

%late
figure(auc_fig_press_late);
subplot(nbins,3,2+((iAuc-1)*3))
hVvsInv = cdfplot(resp_v_late);
hVvsInv.Color = 'g';
hold on
hInvvsV = cdfplot(resp_a_late);
hInvvsV.Color = 'k';
[h p] = kstest2(resp_v_late,resp_a_late);
xlim([-0.05 0.05])
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
title({['ROC bin ' num2str(iAuc)];['late; p= ' num2str(p)]})
axis square

%time-courses
%early
baseStimFrames_press = 0:cycTime:size(tc_vis_early,1)-1;
ttMs_press = (-pre_event_frames:minShortFrames-1)/(cycTime/cycTimeMs);
figure(auc_fig_press_early)
subplot(nbins,3,3+((iAuc-1)*3))
tc1 = mean(tc_vis_early,2)-mean(mean(tc_vis_early(pre_win,:),2));
err1 = std(tc_vis_early,[],2)/sqrt(size(tc_vis_early,2));
shadedErrorBar(ttMs_press,tc1,err1,'g');
hold on
tc2 = mean(tc_aud_early,2)-mean(mean(tc_aud_early(pre_win,:),2));
err2 = std(tc_aud_early,[],2)/sqrt(size(tc_aud_early,2));
shadedErrorBar(ttMs_press,tc2,err2,'k');
hold on
xlim([-15 minShortFrames-1]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames_press/(cycTime/cycTimeMs),':k')
hold on
vline([early_win(1)-pre_event_frames early_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
title({['ROC bin ' num2str(iAuc)]});
% axis square

%late
baseStimFrames_press = 0:cycTime:size(tc_vis_late,1)-1;
ttMs_press = (-pre_event_frames:minLongFrames-1)/(cycTime/cycTimeMs);
figure(auc_fig_press_late)
subplot(nbins,3,3+((iAuc-1)*3))
tc1 = mean(tc_vis_late,2)-mean(mean(tc_vis_late(pre_win,:),2));
err1 = std(tc_vis_late,[],2)/sqrt(size(tc_vis_late,2));
shadedErrorBar(ttMs_press,tc1,err1,'g');
hold on
tc2 = mean(tc_aud_late,2)-mean(mean(tc_aud_late(pre_win,:),2));
err2 = std(tc_aud_late,[],2)/sqrt(size(tc_aud_late,2));
shadedErrorBar(ttMs_press,tc2,err2,'k');
hold on
xlim([-15 minLongFrames-1]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames_press/(cycTime/cycTimeMs),':k')
hold on
vline([late_win(1)-pre_event_frames late_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
title({['ROC bin ' num2str(iAuc)]});
% axis square

end
figure(valOverlay)
legend(legData, legendInfo)
print([fnout '_tc_auROCbins_valOverlay'],'-dpdf','-fillpage')

%fig val and inv tc aligned to last base stim, all cells
tc1 = mean(tc_val_lastBaseAlign,2)-mean(mean(tc_val_lastBaseAlign(pre_win-cycTime,:),1),2);
err1 = std(tc_val_lastBaseAlign,[],2)/sqrt(size(tc_val_lastBaseAlign,2));
tc2 = mean(tc_inv_lastBaseAlign,2)-mean(mean(tc_inv_lastBaseAlign(pre_win-cycTime,:),1),2);
err2 = std(tc_inv_lastBaseAlign,[],2)/sqrt(size(tc_inv_lastBaseAlign,2));

lastBaseAlignFig = figure;
clear legData
l = shadedErrorBar(ttMs,tc1,err1);
l.mainLine.Color = 'k';
l.mainLine.LineWidth = 3;
legData(1) = l.mainLine;
hold on
hold on
l = shadedErrorBar(ttMs,tc2,err2);
l.mainLine.Color = 'c';
l.mainLine.LineWidth = 3;
legData(2) = l.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.02])
xlabel('time (ms)')
ylabel('dF/F')
title('all cells, tc aligned to last base stim')
legend(legData,{'valid';'invalid'});
print([fnout '_tc_lastBaseAlign_val-inv-all'],'-dpdf','-fillpage')


figure(auc_fig);
print([fnout '_auROCBins_val-inv-all'],'-dpdf','-fillpage')
figure(auc_fig_press_late);
print([fnout '_auROCBins_late'],'-dpdf','-fillpage')
figure(auc_fig_press_early);
print([fnout '_auROCBins_early'],'-dpdf','-fillpage')
figure(heatmap_fig);
print([fnout '_heatmap_auROCBins_val-inv-all'],'-dpdf','-fillpage')

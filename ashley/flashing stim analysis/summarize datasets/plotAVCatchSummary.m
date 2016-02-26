function plotAVCatchSummary(datasetStr,cellsInd)
% takes mouse FSAV Ca structure and plots responses to anticipation phase
% of task
close all
av = behavParamsAV;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
if strcmp(rc.name,'ashley')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
cycTime = mouse(1).expt(1).info(1).cyc_time;
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name]; 
iav = 1; % for visual trials
% parameters for accessing data in mouse struct
catchAlign = 3;
targetAlign = 2;
cFA = 3;
cCR = 4;
hits = 1;
misses = 2;

% cycTime = mouse(1).expt(1).info(1).cyc_time;

fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' date '_' mouse_str]); %% maybe lose mouse_str

nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
n = ceil(sqrt(nexp+1));
if (n^2)-n > nexp+1
    n2 = n-1;
else
    n2 = n;
end

%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% set(0,'DefaultaxesFontSize', 16)

tt = -pre_event_frames:post_event_frames-1;
baseStimFrames = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% avg resp accross all direction for success,miss,catch-FA,catch-CR - target responsive cells

respTvsCFig = figure;
i = 1;
resp_tar_hits = [];
resp_tar_misses = [];
resp_catch_fa = [];
resp_catch_cr = [];
% resp_FA = [];
% resp_CR = [];
n_cFA = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_cCR = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_hits = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_misses = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(respTvsCFig);
        subplot(n,n2,i)

        respTrace = [];
        if mouse(imouse).expt(iexp).info.isCatch
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            respTrace(1) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(:,cell_ind,:),2),3), 'c');
            n_cFA(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(:,cell_ind,:),3);
            hold on
            respTrace(2) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cCR).resp(:,cell_ind,:),2),3), 'b');
            n_cCR(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cCR).resp(:,cell_ind,:),3);
            hold on
            respTrace(3) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(:,cell_ind,:),2),3), 'k');
            n_hits(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(:,cell_ind,:),3);
            hold on
            respTrace(4) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(:,cell_ind,:),2),3), 'r');
            n_misses(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(:,cell_ind,:),3);
            hold on
            vline(baseStimFrames,':k')
            hold on
            vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
            hold on
            vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
            hold on
            xlim([-10 20])
            legend(respTrace,{[num2str(n_cFA(i)) ' cFA'];[num2str(n_cCR(i)) ' cCR'];[num2str(n_hits(i)) ' hits'];[num2str(n_misses(i)) ' misses']},'Location','northwest');
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
        
            resp_tar_hits = cat(2, resp_tar_hits, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(:,cell_ind,:),3));
            resp_tar_misses  = cat(2, resp_tar_misses, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(:,cell_ind,:),3));
            resp_catch_fa = cat(2, resp_catch_fa, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(:,cell_ind,:),3));
            resp_catch_cr = cat(2, resp_catch_cr, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cCR).resp(:,cell_ind,:),3));
        end
         i = i+1;       
    end
end

n_cFA_all = sum(n_cFA);
n_cCR_all = sum(n_cCR);
n_misses_all = sum(n_misses);
n_hits_all = sum(n_hits);

respTrace = [];
subplot(n,n2,i)
respTrace(3) = plot(tt, nanmean(resp_tar_hits,2), 'k');
hold on
respTrace(4) = plot(tt, nanmean(resp_tar_misses,2), 'r');
hold on
respTrace(1) = plot(tt, nanmean(resp_catch_fa,2), 'c');
hold on
respTrace(2) = plot(tt, nanmean(resp_catch_cr,2), 'b');
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 20])
legend(respTrace,{[num2str(n_cFA_all) ' cFA'];[num2str(n_cCR_all) ' cCR'];[num2str(n_hits_all) ' hits'];[num2str(n_misses_all) ' misses']},'Location','northwest');
title(['All cells; n = ' num2str(size(resp_tar_hits,2))])
suptitle(titleStr)
figure(respTvsCFig);
print([fnout 'catch_align_TCs' datasetStr '.pdf'], '-dpdf')

%% scatter transient response C-FA vs. T-suc; C-CR vs. T-miss;C-CR vs T-suc; T-miss vs T-suc
scatFAvsH = figure;
scatCRvsM = figure;
scatCRvsH = figure;
scatMvsH = figure;
scatFAvsCR = figure;
i = 1;
resp_h_all = [];
resp_m_all = [];
resp_fa_all = [];
resp_cr_all = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            resp_h = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(pre_win,cell_ind,:),3),1));
            resp_m = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(pre_win,cell_ind,:),3),1));
            resp_fa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(pre_win,cell_ind,:),3),1));
            resp_cr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cCR).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(cFA).resp(pre_win,cell_ind,:),3),1));
            
            resp_h_all = cat(2,resp_h_all,resp_h);
            resp_m_all = cat(2,resp_m_all,resp_m);
            resp_fa_all = cat(2,resp_fa_all,resp_fa);
            resp_cr_all = cat(2,resp_cr_all,resp_cr);
            
            figure(scatFAvsH);
            subplot(n,n2,i)
            scatter(resp_h,resp_fa,50,'k.')
            hold on
            errorbarxy(mean(resp_h),mean(resp_fa),std(resp_h)/length(resp_h),std(resp_fa)/length(resp_fa),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('hits')
            ylabel('fa')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits(i)) ' hits;' num2str(n_cFA(i)) ' cFA']})
            
            figure(scatCRvsM);
            subplot(n,n2,i)
            scatter(resp_m,resp_cr,50,'k.')
            hold on
            errorbarxy(mean(resp_m),mean(resp_cr),std(resp_m)/length(resp_m),std(resp_cr)/length(resp_cr),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('misses')
            ylabel('cr')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_misses(i)) ' misses;' num2str(n_cCR(i)) ' cCR']})
                        
            figure(scatCRvsH);
            subplot(n,n2,i)
            scatter(resp_h,resp_cr,50,'k.')
            hold on
            errorbarxy(mean(resp_h),mean(resp_cr),std(resp_h)/length(resp_h),std(resp_cr)/length(resp_cr),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('hits')
            ylabel('cr')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits(i)) ' hits;' num2str(n_cCR(i)) ' cCR']})
            
            figure(scatMvsH);
            subplot(n,n2,i)
            scatter(resp_h,resp_m,50,'k.')
            hold on
            errorbarxy(mean(resp_h),mean(resp_m),std(resp_h)/length(resp_h),std(resp_m)/length(resp_m),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('hits')
            ylabel('miss')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits(i)) ' hits;' num2str(n_misses(i)) ' misses']})
            
            figure(scatFAvsCR);
            subplot(n,n2,i)
            scatter(resp_cr,resp_fa,50,'k.')
            hold on
            errorbarxy(mean(resp_cr),mean(resp_fa),std(resp_cr)/length(resp_cr),std(resp_fa)/length(resp_fa),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('cr')
            ylabel('fa')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_cCR(i)) ' cCR;' num2str(n_cFA(i)) ' cFA']})
        end
        i = i+1;
    end
end
        
figure(scatFAvsH);
subplot(n,n2,i)
scatter(resp_h_all,resp_fa_all,50,'k.')
hold on
errorbarxy(mean(resp_h_all),mean(resp_fa_all),std(resp_h_all)/length(resp_h_all),std(resp_fa_all)/length(resp_fa_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('hits')
ylabel('fa')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits_all) ' hits;' num2str(n_cFA_all) ' cFA']})
suptitle(titleStr)

figure(scatCRvsM);
subplot(n,n2,i)
scatter(resp_m_all,resp_cr_all,50,'k.')
hold on
errorbarxy(mean(resp_m_all),mean(resp_cr_all),std(resp_m_all)/length(resp_m_all),std(resp_cr_all)/length(resp_cr_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('misses')
ylabel('cr')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(n_misses_all) ' misses;' num2str(n_cCR_all) ' cCR']})
suptitle(titleStr)

figure(scatCRvsH);
subplot(n,n2,i)
scatter(resp_h_all,resp_cr_all,50,'k.')
hold on
errorbarxy(mean(resp_h_all),mean(resp_cr_all),std(resp_h_all)/length(resp_h_all),std(resp_cr_all)/length(resp_cr_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('hits')
ylabel('cr')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits_all) ' hits;' num2str(n_cCR_all) ' cCR']})
suptitle(titleStr)

figure(scatMvsH);
subplot(n,n2,i)
scatter(resp_h_all,resp_m_all,50,'k.')
hold on
errorbarxy(mean(resp_h_all),mean(resp_m_all),std(resp_h_all)/length(resp_h_all),std(resp_m_all)/length(resp_m_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('hits')
ylabel('miss')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hits_all) ' hits;' num2str(n_misses_all) ' misses']})
suptitle(titleStr)

figure(scatFAvsCR);
subplot(n,n2,i)
scatter(resp_cr_all,resp_fa_all,50,'k.')
hold on
errorbarxy(mean(resp_cr_all),mean(resp_fa_all),std(resp_cr_all)/length(resp_cr_all),std(resp_fa_all)/length(resp_fa_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('cr')
ylabel('fa')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(n_cCR_all) ' cCR;' num2str(n_cFA_all) ' cFA']})
suptitle(titleStr)


figure(scatFAvsH);
print([fnout 'catch_align_FAvsH' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsM);
print([fnout 'catch_align_CRvsM' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsH);
print([fnout 'catch_align_CRvsH' datasetStr '.pdf'], '-dpdf')
figure(scatMvsH);
print([fnout 'catch_align_MvsH' datasetStr '.pdf'], '-dpdf')
figure(scatFAvsCR);
print([fnout 'catch_align_FAvsCR' datasetStr '.pdf'], '-dpdf')

%% plot average 
%  plot avg target trace for random subset of cells, all directions
Dirs = mouse(end).info.allDirs; 
colorsT = brewermap(length(Dirs)+15,'YlGn');
colorindT = [3:2:length(Dirs)+12];
colorsT = colorsT(colorindT(1:length(Dirs)),:);


for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        figure;
        suptitle({mouse(imouse).expt(iexp).date, mouse(imouse).expt(iexp).mouse_name})
        cellTarResp = [];
        dirLegend = [];
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        if length(cell_ind) > 9
            cell_ind = sort(randsample(cell_ind,9));
        end
        for icell = 1:length(cell_ind)
            subplot(3,3,icell)
            exptDirs = mouse(imouse).expt(iexp).visTargets;
            for idir = 1:length(exptDirs)  
                dirColInd = find(Dirs == exptDirs(idir));
                cellTarResp(idir) = plot(tt,mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{idir}(:,cell_ind(icell),:),3),'color',colorsT(dirColInd,:));
                dirLegend{idir} = num2str(exptDirs(idir));
                hold on                
            end
            legend(cellTarResp,dirLegend,'Location','northwest')
            vline(baseStimFrames,':k')
            hold on
            vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
            hold on
            vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
            hold on
            xlim([-10 10])
            ylim([-0.1 0.1])
            title(['cell# ' num2str(cell_ind(icell))])
        end
    end
end

end


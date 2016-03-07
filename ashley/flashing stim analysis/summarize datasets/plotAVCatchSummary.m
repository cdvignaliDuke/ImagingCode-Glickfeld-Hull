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
hitmatch = 5;
missmatch = 6;
famatch = 7;
crmatch = 10;

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
suptitle(titleStr)
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
            % invalid hits matched to invalid miss
            respTrace(1) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),2),3), 'c');
            n_cFA(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),3);
            hold on
            %invalid miss matched to invalid hits
            respTrace(2) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),2),3), 'b');
            n_cCR(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),3);
            hold on
            %valid hits matched to valid miss
            respTrace(3) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),2),3), 'k');
            n_hits(i) = size(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),3);
            hold on
            %valid miss match to valid hits
            respTrace(4) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),2),3), 'r');
            n_misses(i) = size(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),3);
            hold on
            vline(baseStimFrames,':k')
            hold on
            vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
            hold on
            vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
            hold on
            xlim([-10 20])
            legend(respTrace,{[num2str(n_cFA(i)) ' invalid hit'];[num2str(n_cCR(i)) ' invalid miss'];[num2str(n_hits(i)) ' hit'];[num2str(n_misses(i)) ' miss']},'Location','northwest');
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
        
            resp_tar_hits = cat(2, resp_tar_hits, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),3));
            resp_tar_misses  = cat(2, resp_tar_misses, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),3));
            resp_catch_fa = cat(2, resp_catch_fa, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),3));
            resp_catch_cr = cat(2, resp_catch_cr, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),3));
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
figure(respTvsCFig);
print([fnout 'catch_align_TCs' datasetStr '.pdf'], '-dpdf')

%% scatter transient response C-FA vs. T-suc; C-CR vs. T-miss;C-CR vs T-suc; T-miss vs T-suc
scatFAvsH = figure; %oucomes 1h 5fa
scatCRvsM = figure;%oucomes 4m 9cr
scatCRvsH = figure;%oucomes 2h 8cr
scatMvsH = figure; % align = 2; oucomes 5h 6m
scatCRvsFA = figure; %oucomes 7fa 10cr
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
n_hvsfa = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_favsh = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_crvsm = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_crvsh = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_mvsh = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
n_crvsfa = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));


for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            
            % fa vs hits matched
            hits = 1;
            fas = 5;
            resp_hvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
            resp_favsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
            resp_hvsfa_all = cat(2,resp_hvsfa_all,resp_hvsfa);
            resp_favsh_all = cat(2,resp_favsh_all,resp_favsh);
            n_hvsfa(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3);
            n_favsH(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3);
            figure(scatFAvsH);
            subplot(n,n2,i)
            scatter(resp_hvsfa,resp_favsh,50,'k.')
            hold on
            errorbarxy(mean(resp_hvsfa),mean(resp_favsh),std(resp_hvsfa)/length(resp_hvsfa),std(resp_favsh)/length(resp_favsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('invalid hits (fa)')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hvsfa(i)) ' valid hits;' num2str(n_favsH(i)) ' invalid hits']})
            
            % cr vs miss matched
            misses = 4;
            crs = 9;
            resp_mvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
            resp_crvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_mvscr_all = cat(2,resp_mvscr_all,resp_mvscr);
            resp_crvsm_all = cat(2,resp_crvsm_all,resp_crvsm);
            n_crvsm(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3);
            figure(scatCRvsM);
            subplot(n,n2,i)
            scatter(resp_mvscr,resp_crvsm,50,'k.')
            hold on
            errorbarxy(mean(resp_mvscr),mean(resp_crvsm),std(resp_mvscr)/length(resp_mvscr),std(resp_crvsm)/length(resp_crvsm),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid miss')
            ylabel('invalid miss')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsm(i)) ' valid miss;' num2str(n_crvsm(i)) ' invalid miss']})
                  
            % cr vs hits matched
            hits = 2;
            crs = 8;
            resp_hvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
            resp_crvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_hvscr_all = cat(2,resp_hvscr_all,resp_hvscr);
            resp_crvsh_all = cat(2,resp_crvsh_all,resp_crvsh);
            n_crvsh(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3);
            figure(scatCRvsH);
            subplot(n,n2,i)
            scatter(resp_hvscr,resp_crvsh,50,'k.')
            hold on
            errorbarxy(mean(resp_hvscr),mean(resp_crvsh),std(resp_hvscr)/length(resp_hvscr),std(resp_crvsh)/length(resp_crvsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('invalid miss')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsh(i)) ' valid hits;' num2str(n_crvsh(i)) ' invalid miss']})
        
            % hits vs miss matched
            hits = 5;
            misses = 6;
            resp_hvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
            resp_mvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
            resp_hvsm_all = cat(2,resp_hvsm_all,resp_hvsm);
            resp_mvsh_all = cat(2,resp_mvsh_all,resp_mvsh);
            n_mvsh(i) = size(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3);
            figure(scatMvsH);
            subplot(n,n2,i)
            scatter(resp_hvsm,resp_mvsh,50,'k.')
            hold on
            errorbarxy(mean(resp_hvsm),mean(resp_mvsh),std(resp_hvsm)/length(resp_hvsm),std(resp_mvsh)/length(resp_mvsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('valid miss')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_mvsh(i)) ' valid hits;' num2str(n_mvsh(i)) ' valid miss']})
            
            % fa vs cr matched
            fas = 7;
            crs = 10;
            resp_favscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
            resp_crvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_favscr_all = cat(2,resp_favscr_all,resp_favscr);
            resp_crvsfa_all = cat(2,resp_crvsfa_all,resp_crvsfa);
            n_crvsfa(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3);
            figure(scatCRvsFA);
            subplot(n,n2,i)
            scatter(resp_favscr,resp_crvsfa,50,'k.')
            hold on
            errorbarxy(mean(resp_favscr),mean(resp_crvsfa),std(resp_favscr)/length(resp_favscr),std(resp_crvsfa)/length(resp_crvsfa),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('invalid hits')
            ylabel('invalid miss')
            title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsfa(i)) ' invalid hits;' num2str(n_crvsfa(i)) ' invalid miss']})
        end
        i = i+1;
    end
end
        
figure(scatFAvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvsfa_all,resp_favsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvsfa_all),mean(resp_favsh_all),std(resp_hvsfa_all)/length(resp_hvsfa_all),std(resp_favsh_all)/length(resp_favsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('invalid hits')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_hvsfa)) ' valid hits;' num2str(sum(n_favsh)) ' invalid hits']})

figure(scatCRvsM);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_mvscr_all,resp_crvsm_all,50,'k.')
hold on
errorbarxy(mean(resp_mvscr_all),mean(resp_crvsm_all),std(resp_mvscr_all)/length(resp_mvscr_all),std(resp_crvsm_all)/length(resp_crvsm_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid miss')
ylabel('invalid miss')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsm)) ' valid miss;' num2str(sum(n_crvsm)) ' invalid miss']})

figure(scatCRvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvscr_all,resp_crvsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvscr_all),mean(resp_crvsh_all),std(resp_hvscr_all)/length(resp_hvscr_all),std(resp_crvsh_all)/length(resp_crvsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('invalid miss')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsh)) 'valid hits;' num2str(sum(n_crvsh)) ' invalid miss']})

figure(scatMvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvsm_all,resp_mvsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvsm_all),mean(resp_mvsh_all),std(resp_hvsm_all)/length(resp_hvsm_all),std(resp_mvsh_all)/length(resp_mvsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('valid miss')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_mvsh)) ' valid hits;' num2str(sum(n_mvsh)) ' valid miss']})

figure(scatCRvsFA);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_favscr_all,resp_crvsfa_all,50,'k.')
hold on
errorbarxy(mean(resp_favscr_all),mean(resp_crvsfa_all),std(resp_favscr_all)/length(resp_favscr_all),std(resp_crvsfa_all)/length(resp_crvsfa_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('invalid hit')
ylabel('invalid miss')
title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsfa)) ' invalid hit;' num2str(sum(n_crvsfa)) ' invalid miss']})


figure(scatFAvsH);
print([fnout 'catch_align_FAvsH' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsM);
print([fnout 'catch_align_CRvsM' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsH);
print([fnout 'catch_align_CRvsH' datasetStr '.pdf'], '-dpdf')
figure(scatMvsH);
print([fnout 'catch_align_MvsH' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsFA);
print([fnout 'catch_align_CRvsFA' datasetStr '.pdf'], '-dpdf')

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


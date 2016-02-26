function plotAVAnticipationSummary(datasetStr,cellsInd)
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
minTrialLengthFrames = mouse(1).expt(1).info.minTrialLengthFrames;
ialign = 1;
cycTime = mouse(1).expt(1).info(1).cyc_time;
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name]; 

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

tt = -pre_event_frames:minTrialLengthFrames-1;
baseStimFrames = 0:cycTime:minTrialLengthFrames-1;

%% plot avg trace for each expt and avg across all expt - success and responsive cells only
respAvsVFig = figure;
i = 1;
resp_vis = [];
resp_aud = [];
% resp_FA = [];
% resp_CR = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(respAvsVFig);
        subplot(n,n2,i)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),2),3), 'g');
        hold on
        plot(tt,mean(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),2),3), 'k');
        hold on
%         vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%         vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        vline(baseStimFrames,':k')
        xlim([-10 minTrialLengthFrames])
        ylim([-0.01 0.03]);
        resp_vis = cat(2, resp_vis, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),3));
        resp_aud = cat(2, resp_aud, mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),3));
%         resp_FA = cat(2, resp_FA, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(3).cmlvResp(:,cell_ind,:),3));
%         resp_CR = cat(2, resp_CR, mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(4).cmlvResp(:,cell_ind,:),3));
        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
    end
end
subplot(n,n2,i)
plot(tt, nanmean(resp_vis,2), 'g')
hold on
plot(tt, nanmean(resp_aud,2), 'k')
% plot(tt, nanmean(resp_FA,2), 'c')
% plot(tt, nanmean(resp_CR,2), 'b')
% vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
% vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
xlim([-10 minTrialLengthFrames])
ylim([-0.01 0.03]);
vline(baseStimFrames,':k')
ylim([-0.01 0.05]);
title(['All cells; n = ' num2str(size(resp_vis,2))])
suptitle({titleStr, 'Hits: Black; FAs: Cyan; CR: Blue; Auditory: Green'})
figure(respAvsVFig);
print([fnout 'press_align_TCs' datasetStr '.pdf'], '-dpdf')


%% plot avg trace for random subset of cells

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        figure;
        suptitle({mouse(imouse).expt(iexp).date, mouse(imouse).expt(iexp).mouse_name})
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        if length(cell_ind) > 9
            cell_ind = sort(randsample(cell_ind,9));
        end
        for icell = 1:length(cell_ind)
            subplot(3,3,icell)
            plot(tt,mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind(icell),:),3),'g')
            hold on
            plot(tt,mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind(icell),:),3),'k')
            hold on
            vline(baseStimFrames,':k')
            title(['cell# ' num2str(cell_ind(icell))])
        end
    end
end

%% plot scatter of integral for each experiment, across experiments;
resp_vis = [];
resp_aud = [];
intAvsVFig = figure;
i = 1;

for imouse = 1:size(mouse,2)
    resp_vis_mouse{imouse} = [];
    resp_aud_mouse{imouse} = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(intAvsVFig);
        subplot(n,n2,i)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        % add STE color axis here
        scatter(iV,iA,50,'k.');
        hold on
        errorbarxy(mean(iV),mean(iA),std(iV)/length(iV),std(iA)/length(iA),{'ro','r','r'});
%         xlim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
%         ylim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
            xlim([-5 15])
            ylim([-5 15])
%         
        hold on
        plot([-10:0.1:20],[-10:0.1:20],'k--')
        axis square
        xlabel('Vis Tr Resp')
        ylabel('Aud Tr Resp')
        
        resp_vis = cat(2, resp_vis, iV);
        resp_aud = cat(2, resp_aud, iA);
        
        resp_vis_mouse{imouse} = cat(2,resp_vis_mouse{imouse}, iV);
        resp_aud_mouse{imouse} = cat(2,resp_aud_mouse{imouse}, iA);

        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
    end
end

subplot(n,n2,i)
scatter(resp_vis,resp_aud,50,'k.')
hold on
errorbarxy(mean(resp_vis),mean(resp_aud),std(resp_vis)/length(resp_vis),std(resp_aud)/length(resp_aud),{'ro','r','r'})
hold on
for imouse = 1:size(mouse,2)
   errorbarxy(mean(resp_vis_mouse{imouse}),mean(resp_aud_mouse{imouse}),std(resp_vis_mouse{imouse})/length(resp_vis_mouse{imouse}),std(resp_aud_mouse{imouse})/length(resp_vis_mouse{imouse}), {[av(mouse_ind(imouse)).col_str 'o'],av(mouse_ind(imouse)).col_str,av(mouse_ind(imouse)).col_str});
   hold on
   mouse_mean(imouse) = scatter(mean(resp_vis_mouse{imouse}),mean(resp_aud_mouse{imouse}),av(mouse_ind(imouse)).col_str,'o', 'filled');
   %    mouseLegend{imouse} = str{imouse};
   hold on
end
legend(mouse_mean,str,'Location','southeast')
% xlim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
% ylim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
xlim([-5 15])
ylim([-5 15])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

title(['All cells; n = ' num2str(size(resp_vis,2))])
suptitle(titleStr)
figure(intAvsVFig);
print([fnout 'press_align_int' datasetStr '.pdf'], '-dpdf')

%% plot scatter of integral for each experiment, across experiments; ORI SELECT LABELLED 
resp_vis = cell([1,4]);
resp_aud = cell([1,4]);
intAvsVFig_ori = figure;
i = 1;
oriColMat = ['k';'g';'b';'r'];

for imouse = 1:size(mouse,2)
%     resp_vis_mouse{imouse} = [];
%     resp_aud_mouse{imouse} = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        figure(intAvsVFig_ori);
        subplot(n,n2,i)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        for iori = 1:4
            ori_ind = intersect(mouse(imouse).expt(iexp).cells(iori+1).ind,cell_ind);
            iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,ori_ind,:)),3));
            iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,ori_ind,:)),3));
            % add STE color axis here
            scatter(iV,iA,50,[oriColMat(iori) '.']);
            hold on
%             errorbarxy(mean(iV),mean(iA),std(iV)/length(iV),std(iA)/length(iA),{[oriColMat(iori) 'o'],oriColMat(iori),oriColMat(iori)});
            oriMean(iori) = scatter(mean(iV),mean(iA),oriColMat(iori),'o','filled');
            if length([iV iA]) > 0;
%             xlim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
%             ylim([floor(min([iV iA]))-(floor(min([iV iA]))/2) ceil(max([iV iA]))+(ceil(max([iV iA]))/2)]);
            xlim([-5 15])
            ylim([-5 15])
            end
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('Vis Tr Resp')
            ylabel('Aud Tr Resp')
            
            
            resp_vis{iori} = cat(2, resp_vis{iori}, iV);
            resp_aud{iori} = cat(2, resp_aud{iori}, iA);
        end

        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
    end
end

subplot(n,n2,i)
% hold on
% errorbarxy(mean(resp_vis),mean(resp_aud),std(resp_vis)/length(resp_vis),std(resp_aud)/length(resp_aud),{'ro','r','r'})
% hold on
oriLegend = {'0','45','90','135'};
for iori = 1:4
    scatter(resp_vis{iori},resp_aud{iori},50,oriColMat(iori),'.')
    hold on
    errorbarxy(mean(resp_vis{iori}),mean(resp_aud{iori}),std(resp_vis{iori})/length(resp_vis{iori}),std(resp_aud{iori})/length(resp_aud{iori}), {[oriColMat(iori) 'o'],oriColMat(iori),oriColMat(iori)});
    hold on
    ori_mean(iori) = scatter(mean(resp_vis{iori}),mean(resp_aud{iori}),oriColMat(iori),'o', 'filled');
    hold on
end
legend(ori_mean,oriLegend,'Location','southeast')
% xlim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
% ylim([floor(min([resp_vis resp_aud]))-(floor(min([resp_vis resp_aud]))/2) ceil(max([resp_vis resp_aud]))+(ceil(max([resp_vis resp_aud]))/2)]);
xlim([-5 15])
ylim([-5 15])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

title(['All cells; n = ' num2str(size(resp_vis,2))])
suptitle(titleStr)
figure(intAvsVFig_ori);
print([fnout 'press_align_intOri' datasetStr '.pdf'], '-dpdf')

%% plot histogram of number of excited, inhibited, and non-responsive cells by resp integral

histEvsIvsNRcells_V = figure;
histEvsIvsNRcells_A = figure;
nbins = 100;
i = 1;
resp_vis = [];
resp_aud = [];
exc_ind_all = [];
inh_ind_all = [];
% resp_FA = [];
% resp_CR = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        
        cell_ind = 1:mouse(imouse).expt(iexp).info.nCells;
        exc_ind = mouse(imouse).expt(iexp).cells(6).ind;
        inh_ind = mouse(imouse).expt(iexp).cells(7).ind;
        
        iV = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        iA = squeeze(mean(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:end,cell_ind,:)),3));
        
        bin_edges = linspace(min([iV iA]),max([iV iA]),nbins);
        
        [barVexc barVexcInd] = histc(iV(exc_ind),bin_edges);
        [barVinh barVinhInd] = histc(iV(inh_ind),bin_edges);
        [barVnr barVnrInd] = histc(iV(setdiff(cell_ind, [exc_ind' inh_ind'])),bin_edges);
        
        [barAexc barAexcInd] = histc(iA(exc_ind),bin_edges);
        [barAinh barAinhInd] = histc(iA(inh_ind),bin_edges);
        [barAnr barAnrInd] = histc(iA(setdiff(cell_ind, [exc_ind' inh_ind'])),bin_edges);
        
        figure(histEvsIvsNRcells_V);
        subplot(n,n2,i)
        hE_V = bar(bin_edges,barVexc,'r');
        hE_V_ch = get(hE_V,'child');
        set(hE_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        hI_V = bar(bin_edges,barVinh,'b');
        hI_V_ch = get(hI_V,'child');
        set(hI_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        hNR_V = bar(bin_edges,barVnr,'k');
        hNR_V_ch = get(hNR_V,'child');
        set(hNR_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        vline(0,'k')
        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
        xlabel('resp int')
        ylabel('# cells')
        
        figure(histEvsIvsNRcells_A);
        subplot(n,n2,i)
        hE_A = bar(bin_edges,barAexc,'r');
        hE_A_ch = get(hE_A,'child');
        set(hE_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        hI_A = bar(bin_edges,barVinh,'b');
        hI_A_ch = get(hI_A,'child');
        set(hI_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        hNR_A = bar(bin_edges,barVnr,'k');
        hNR_A_ch = get(hNR_A,'child');
        set(hNR_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
        hold on
        vline(0,'k')
        title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
        xlabel('resp int')
        ylabel('# cells')
        
        resp_vis = cat(2, resp_vis, iV);
        resp_aud = cat(2, resp_aud, iA);
        
        if i > 1
            exc_ind_all = cat(1,exc_ind_all,exc_ind+runCells);
            inh_ind_all = cat(1,inh_ind_all,inh_ind+runCells);
            runCells = length(cell_ind)+runCells;
        else
            exc_ind_all = exc_ind;
            inh_ind_all = inh_ind;
            runCells = length(cell_ind);
        end
        
        i = i+1;
    end
end

bin_edges = linspace(min([resp_vis resp_aud]),max([resp_vis resp_aud]),nbins);
        
[barVexc barVexcInd] = histc(resp_vis(exc_ind_all),bin_edges);
[barVinh barVinhInd] = histc(resp_vis(inh_ind_all),bin_edges);
[barVnr barVnrInd] = histc(resp_vis(setdiff(1:runCells, [exc_ind_all' inh_ind_all'])),bin_edges);

[barAexc barAexcInd] = histc(resp_aud(exc_ind_all),bin_edges);
[barAinh barAinhInd] = histc(resp_aud(inh_ind_all),bin_edges);
[barAnr barAnrInd] = histc(resp_aud(setdiff(1:runCells, [exc_ind_all' inh_ind_all'])),bin_edges);

figure(histEvsIvsNRcells_V);
subplot(n,n2,i)
hE_V = bar(bin_edges,barVexc,'r');
hE_V_ch = get(hE_V,'child');
set(hE_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
hold on
hI_V = bar(bin_edges,barVinh,'b');
hI_V_ch = get(hI_V,'child');
set(hI_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
hold on
hNR_V = bar(bin_edges,barVnr,'k');
hNR_V_ch = get(hNR_V,'child');
set(hNR_V_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         set(gca,'XTick',bin_edges)
hold on
vline(0,'k')
title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
xlabel('resp int')
ylabel('# cells')

figure(histEvsIvsNRcells_A);
subplot(n,n2,i)
hE_A = bar(bin_edges,barAexc,'r');
hE_A_ch = get(hE_A,'child');
set(hE_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
hold on
hI_A = bar(bin_edges,barVinh,'b');
hI_A_ch = get(hI_A,'child');
set(hI_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
hold on
hNR_A = bar(bin_edges,barVnr,'k');
hNR_A_ch = get(hNR_A,'child');
set(hNR_A_ch,'FaceAlpha',0.25,'EdgeAlpha',0);
%         set(gca,'XTick',bin_edges)
hold on
vline(0,'k')
title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
xlabel('resp int')
ylabel('# cells')

figure(histEvsIvsNRcells_V);
print([fnout 'press_align_Vint_hist' datasetStr '.pdf'], '-dpdf')
figure(histEvsIvsNRcells_A);
print([fnout 'press_align_Aint_hist' datasetStr '.pdf'], '-dpdf')
end

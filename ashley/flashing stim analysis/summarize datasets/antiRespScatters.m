alltr = 3;

resp.av(visual).all = [];
resp.av(auditory).all = [];
resp.av(visual).early = [];
resp.av(auditory).early = [];
resp.av(visual).late = [];
resp.av(auditory).late = [];
resp.ori(1).name = '0';
resp.ori(2).name = '45';
resp.ori(3).name = '90';
resp.ori(4).name = '135';
resp.ori(5).name = 'non-slctv';
resp.ori(1).ind = [];
resp.ori(2).ind = [];
resp.ori(3).ind = [];
resp.ori(4).ind = [];
resp.ori(5).ind = [];
resp.av(alltr).tc = [];
resp.av(alltr).all = [];
resp.av(alltr).early = [];
resp.av(alltr).late = [];

ncells_ind = 0;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        slctv = [];
        for iori = 1:4
            ori_ind = find(ismember(cell_ind,mouse(imouse).expt(iexp).cells(iori+1).ind));
            resp.ori(iori).ind = cat(1,resp.ori(iori).ind,ori_ind+ncells_ind);
            slctv = cat(1,slctv,ori_ind);
            if iori == 4
                resp.ori(5).ind = cat(2,resp.ori(5).ind,setdiff(1:length(cell_ind),slctv)+ncells_ind);
            end
        end
        v = mean(mouse(imouse).expt(iexp).align(ialign).av(visual).outcome(hits).cmlvCycResp{earlyCycMax},3);
        a = mean(mouse(imouse).expt(iexp).align(ialign).av(auditory).outcome(hits).cmlvCycResp{earlyCycMax},3);
                
        alltr_temp = mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(visual).outcome(hits).cmlvCycResp{earlyCycMax}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(auditory).outcome(hits).cmlvCycResp{earlyCycMax}(:,cell_ind,:)),3);
        resp.av(alltr).early = cat(2,resp.av(alltr).all,trapz(alltr_temp(early_win,:)));
        
        resp.av(visual).early = cat(2,resp.av(visual).early,trapz(v(early_win,cell_ind)));
        resp.av(auditory).early = cat(2,resp.av(auditory).early,trapz(a(early_win,cell_ind)));
        
        v = mean(mouse(imouse).expt(iexp).align(ialign).av(visual).outcome(hits).cmlvCycResp{lateCycMax},3);
        a = mean(mouse(imouse).expt(iexp).align(ialign).av(auditory).outcome(hits).cmlvCycResp{lateCycMax},3);
        alltr_temp = mean(cat(3,mouse(imouse).expt(iexp).align(ialign).av(visual).outcome(hits).cmlvCycResp{lateCycMax}(:,cell_ind,:),mouse(imouse).expt(iexp).align(ialign).av(auditory).outcome(hits).cmlvCycResp{lateCycMax}(:,cell_ind,:)),3);
        
        resp.av(visual).all = cat(2,resp.av(visual).all,trapz(v(pre_event_frames:end,cell_ind)));
        resp.av(auditory).all = cat(2,resp.av(auditory).all,trapz(a(pre_event_frames:end,cell_ind)));
        resp.av(alltr).all = cat(2,resp.av(alltr).all,trapz(alltr_temp(pre_event_frames:end,:)));
        resp.av(visual).late = cat(2,resp.av(visual).late,trapz(v(late_win,cell_ind,:)));
        resp.av(auditory).late = cat(2,resp.av(auditory).late,trapz(a(late_win,cell_ind,:)));
        resp.av(alltr).late = cat(2,resp.av(alltr).early,trapz(alltr_temp(late_win,:)));
        resp.av(alltr).tc = cat(2,resp.av(alltr).tc,alltr_temp);
        
        ncells_ind = length(cell_ind)+ncells_ind;
        
    end
end


%%
cellRespFig = figure;

colors = brewermap(5,'*RdPu');
colors(5,:) = [0.5 0.5 0.5];
ax_range = [-2 5];

subplot(1,3,1)
H = [];
for iori = 1:5
    H = scatter(resp.av(visual).all(:,resp.ori(iori).ind),resp.av(auditory).all(:,resp.ori(iori).ind),15);
    H.MarkerFaceColor = colors(iori,:);
    H.MarkerEdgeColor = 'none';
    hold on    
end  
title('entire trial')
% plot([-5:0.5:5],[-5:0.5:5],'k--')
% xlim(ax_range)
% ylim(ax_range)
% legend({resp.ori.name},'Location','SouthEast')
% axis square

subplot(1,3,2)
H = [];
for iori = 1:5
    H = scatter(resp.av(visual).early(:,resp.ori(iori).ind),resp.av(auditory).early(:,resp.ori(iori).ind),15);
    H.MarkerFaceColor = colors(iori,:);
    H.MarkerEdgeColor = 'none';
    hold on    
end  
title('early phase')
% plot([-5:0.5:5],[-5:0.5:5],'k--')
% xlim(ax_range)
% ylim(ax_range)
% legend({resp.ori.name},'Location','SouthEast')
% axis square

subplot(1,3,3)
H = [];
for iori = 1:5
    H = scatter(resp.av(visual).late(:,resp.ori(iori).ind),resp.av(auditory).late(:,resp.ori(iori).ind),15);
    H.MarkerFaceColor = colors(iori,:);
    H.MarkerEdgeColor = 'none';
    hold on    
end  
title('late phase')

for iplot = 1:3
    subplot(1,3,iplot)
    plot([-5:0.5:5],[-5:0.5:5],'k--')
    xlim(ax_range)
    ylim(ax_range)
    legend({resp.ori.name},'Location','SouthEast')
    axis square
    xlabel('visual')
    ylabel('auditory')
end

print([fnout 'press_align_intOri_scatter' datasetStr '.pdf'], '-dpdf')

%% all trials tc and cdf by ori 

alltr_fig = figure;

for iori = 1:5
   subplot(3,2,iori)
   H = shadedErrorBar(ttMs,mean(resp.av(alltr).tc(1:length(ttMs),resp.ori(iori).ind),2),std(resp.av(alltr).tc(1:length(ttMs),resp.ori(iori).ind),[],2)/sqrt(size(resp.av(alltr).tc(:,resp.ori(iori).ind),2)));
   H.mainLine.Color = colors(iori,:);
   H.mainLine.LineWidth = 1.5;
%    H.patch.FaceColor = H.patch.FaceColor.*colors(iori,:);
   H.edge(1).Color = [1 1 1];
   H.edge(2).Color = [1 1 1];
   ylim([-0.05 .15])
   xlabel('time (ms)')
   ylabel('dF/F')
   title(resp.ori(iori).name)
   vline(baseStimFrames/(cycTime/cycTimeMs),'k--')
end

subplot(3,2,6)
for iori = 1:5
    H = cdfplot(resp.av(alltr).all(resp.ori(iori).ind));
    H.Color = colors(iori,:);
    H.LineWidth = 1.5;
    hold on
end
xlabel('dF/F - int')
ylabel('fraction cells')


print([fnout 'press_align_allTrials_ori_TCandCDF' datasetStr '.pdf'], '-dpdf')
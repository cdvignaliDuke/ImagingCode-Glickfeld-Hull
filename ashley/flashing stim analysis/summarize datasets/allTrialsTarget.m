visual = 1;
auditory = 2;
hits = 1;

dirs = [];
for imouse =1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        dirs = unique(cat(2,dirs,mouse(imouse).expt(iexp).visTargets));
    end
end

cells.tc = cell(1,length(dirs));
cells.ori(1).name = '0';
cells.ori(2).name = '45';
cells.ori(3).name = '90';
cells.ori(4).name = '135';
cells.ori(5).name = 'non-slctv';
cells.ori(1).ind = [];
cells.ori(2).ind = [];
cells.ori(3).ind = [];
cells.ori(4).ind = [];
cells.ori(5).ind = [];
ncells_ind = 0;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(14).ind;
        slctv = [];
        for iori = 1:4
            ori_ind = find(ismember(cell_ind,mouse(imouse).expt(iexp).cells(iori+1).ind));
            cells.ori(iori).ind = cat(1,cells.ori(iori).ind,ori_ind+ncells_ind);
            slctv = cat(1,slctv,ori_ind);
            if iori == 4
                cells.ori(5).ind = cat(2,cells.ori(5).ind,setdiff(1:length(cell_ind),slctv)+ncells_ind);
            end
        end
        
        ncells_ind = length(cell_ind)+ncells_ind;
        
        dirs_temp_ind = find(ismember(dirs,mouse(imouse).expt(iexp).visTargets));
        
        r = cellfun(@(x) mean(x,3),mouse(imouse).expt(iexp).align(ialign).av(visual).outcome(hits).stimResp,'unif',false);
        cells.tc(dirs_temp_ind) = cellfun(@(x,y) cat(2,x,y), cells.tc(dirs_temp_ind), r, 'unif',false);
                
    end
end


%% ori tuned cells

baseStimFrames = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

cellRespFig = figure;

colors = brewermap(length(dirs),'Greens');
colors(1,:) = [0 0 0];

for iori = 1:5
   r = cellfun(@(x) mean(x(:,cells.ori(iori).ind),2)-mean(mean(x(pre_win,cells.ori(iori).ind),2),1),cells.tc,'unif',false);
   r_err = cellfun(@(x) std(x(:,cells.ori(iori).ind),[],2)/sqrt(size(x(:,cells.ori(iori).ind),2)),cells.tc,'unif',false);
   r = cat(2,r{:});
   r_err = cat(2,r_err{:});
   subplot(3,2,iori)
   for i = 1:length(dirs)
       P1 = shadedErrorBar(ttMs,r(:,i),r_err(:,i),'k',1);
       P1.mainLine.Color = colors(i,:);
       P1.mainLine.LineWidth = 1.5;
       P1.edge(1).Color = [1 1 1];
       P1.edge(2).Color = [1 1 1];
       if iori == 1
           P_leg(i) = P1.mainLine;
       end
       hold on
   end
   if iori == 1
       legend(P_leg,strread(num2str(dirs),'%s'),'Location','NorthWest')
   end
   title([cells.ori(iori).name ' slctv'])
%    xlim([-1000 1000])
   ylim([-0.05 0.05])
end
subplot(3,2,6)
r = cellfun(@(x) mean(x,2)-mean(mean(x(pre_win,:),2),1),cells.tc,'unif',false);
r_err = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),cells.tc,'unif',false);
r = cat(2,r{:});
r_err = cat(2,r_err{:});
for i = 1:length(dirs)
       P1 = shadedErrorBar(ttMs,r(:,i),r_err(:,i),'k',1);
       P1.mainLine.Color = colors(i,:);
       P1.mainLine.LineWidth = 1.5;
       P1.edge(1).Color = [1 1 1];
       P1.edge(2).Color = [1 1 1];
       hold on
       P_leg(i) = P1.mainLine;
end
title('all cells')
   xlim([-1000 1000])
ylim([-0.01 0.02])

print([fnout 'target_align_oriTuned_TC' datasetStr '.pdf'], '-dpdf')
%% all cells
figure;
subplot(1,2,1)
r = cellfun(@(x) mean(x,2)-mean(mean(x(pre_win,:),2),1),cells.tc,'unif',false);
r_err = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),cells.tc,'unif',false);
r = cat(2,r{:});
r_err = cat(2,r_err{:});
for i = 1:length(dirs)
       P1 = shadedErrorBar(ttMs,r(:,i),r_err(:,i),'k');
       P1.mainLine.Color = colors(i,:);
       P1.mainLine.LineWidth = 1.5;
       P1.edge(1).Color = [1 1 1];
       P1.edge(2).Color = [1 1 1];
       hold on
       P_leg(i) = P1.mainLine;
end
title('all cells')
ylabel('dF/F')
xlabel('time (ms)')
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.02])
leg_str = cat(1,{'auditory'},strread(num2str(dirs(2:end)),'%s'));
legend(P_leg,leg_str,'Location','NorthWest')
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')

r = cellfun(@(x) mean(x(trans_win,:),1)-mean(x(pre_win,:),1),cells.tc,'unif',false);
subplot(1,2,2)
for i = 1:length(dirs)
    P = cdfplot(r{i})
    P.Color = colors(i,:)
    hold on
end
title('resp to target')
xlim([-0.05 0.05])
xlabel('dF/F')
ylabel('fraction of cells')


print([fnout 'target_align_allcells_TCandCDF' datasetStr '.pdf'], '-dpdf')
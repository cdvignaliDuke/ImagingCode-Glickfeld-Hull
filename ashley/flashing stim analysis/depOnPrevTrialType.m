% must run beginning of plotAVAnticipationSummary 

%% create structure of visual and auditory trial types 
% example: visual trials will be 1x3 cell array that's a cumulative
% response for all trials, prev same trials, and prev different trials,
% corresponding  naming array also made here

vis_mat = cell(1,3);
vis_nTrials = zeros(1,3);
vis_mat_name = cell(1,3);
aud_mat = cell(1,3);
aud_nTrials = zeros(1,3);
aud_mat_name = cell(1,3);
vis_mat_ind = [1 3 4];
aud_mat_ind = [2 5 6];
matfil = 0;
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(respCellsInd).ind,cell_ind);
        
        for imat = 1:3
        vis_mat{imat} = cat(2,vis_mat{imat},mean(mouse(imouse).expt(iexp).align(ialign).av(vis_mat_ind(imat)).outcome(1).cmlvResp(:,cell_ind,:),3));
        aud_mat{imat} = cat(2,aud_mat{imat},mean(mouse(imouse).expt(iexp).align(ialign).av(aud_mat_ind(imat)).outcome(1).cmlvResp(:,cell_ind,:),3));
        
        vis_nTrials(imat) = vis_nTrials(imat)+size(mouse(imouse).expt(iexp).align(ialign).av(vis_mat_ind(imat)).outcome(1).cmlvResp,3);
        aud_nTrials(imat) = aud_nTrials(imat)+size(mouse(imouse).expt(iexp).align(ialign).av(aud_mat_ind(imat)).outcome(1).cmlvResp,3);
        
        if matfil == 0    
        vis_mat_name{imat} = mouse(imouse).expt(iexp).align(ialign).av(vis_mat_ind(imat)).name;
        aud_mat_name{imat} = mouse(imouse).expt(iexp).align(ialign).av(aud_mat_ind(imat)).name;
        end
        end
        matfil = 1;
    end
end

vis_mat_mean = cellfun(@(x) mean(x,2), vis_mat,'unif',false);
aud_mat_mean = cellfun(@(x) mean(x,2), aud_mat,'unif',false);

vis_mat_ste = cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),vis_mat,'unif',false);
aud_mat_ste = cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),aud_mat,'unif',false);

%% plot comparisons

prevTrialAnalysisFig = figure;
suptitle({titleStr; ['cmlv resp by previous trial type']})

% all vis vs all aud
subplot(4,2,1)
l1_ind = 1;
l2_ind = 1;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},aud_mat_name{l1_ind},'Location','northwest')

% vis same vs diff
subplot(4,2,3)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,vis_mat_mean{l2_ind},vis_mat_ste{l2_ind},'g');
l2.mainLine.Color = l2.mainLine.Color./2;
l2.patch.FaceColor = l2.mainLine.Color+0.25;
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
% aud same vs diff
subplot(4,2,4)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,aud_mat_mean{l1_ind},aud_mat_ste{l1_ind},'k');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
l2.mainLine.Color = l2.mainLine.Color+0.5;
l2.patch.FaceColor = l2.patch.FaceColor./2;
hold on
legend([l1.mainLine l2.mainLine],aud_mat_name{l1_ind},aud_mat_name{l2_ind},'Location','northwest')

% vis same vs aud same
subplot(4,2,5)
l1_ind = 2;
l2_ind = 2;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},aud_mat_name{l2_ind},'Location','northwest')
% vis diff vs aud diff
subplot(4,2,6)
l1_ind = 3;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},aud_mat_name{l2_ind},'Location','northwest')

% vis same vs aud diff
subplot(4,2,7)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},aud_mat_name{l2_ind},'Location','northwest')
%vis diff vs aud same
subplot(4,2,8)
l1_ind = 3;
l2_ind = 2;
l1 = shadedErrorBar(ttMs,vis_mat_mean{l1_ind},vis_mat_ste{l1_ind},'g');
hold on
l2 = shadedErrorBar(ttMs,aud_mat_mean{l2_ind},aud_mat_ste{l2_ind},'k');
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},aud_mat_name{l2_ind},'Location','northwest')


figure(prevTrialAnalysisFig)
for iplot = 1:8;
    subplot(4,2,iplot)
    hold on
    xlim([-200 max(ttMs)])
    ylim([-0.01 0.05])
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    xlabel('t(ms)')
    ylabel('dF/F')
end


print([fnout 'press_align_respByPrevTrial' datasetStr '.pdf'], '-dpdf');


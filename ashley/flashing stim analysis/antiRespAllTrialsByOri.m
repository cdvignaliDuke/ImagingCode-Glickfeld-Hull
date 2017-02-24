int_all_ori = cell([1,4]);
resp_all_ori = cell([1,4]);
% intAllFig_ori = figure;
% suptitle(titleStr)
i = 1;
oriColMat = ['k';'g';'b';'r'];
oriStr = {'0';'45';'90';'135'};

for imouse = 1:size(mouse,2)
%     resp_vis_mouse{imouse} = [];
%     resp_aud_mouse{imouse} = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        i = i+1;
        cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
        cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        for iori = 1:4
            ori_ind = intersect(mouse(imouse).expt(iexp).cells(iori+1).ind,cell_ind);
            tempAll = cat(3,mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,ori_ind,:),mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,ori_ind,:));
            iAll = squeeze(mean(trapz(tempAll(pre_event_frames:end,:,:)),3));
            tcAll = mean(tempAll,3); 
            int_all_ori{iori} = cat(2, int_all_ori{iori}, iAll);
            resp_all_ori{iori} = cat(2, resp_all_ori{iori}, tcAll);
        end
    end
end

%cdf plot of tuned cell respones
% figure(intAllFig_ori);
figure(allTrialsFig);
subplot(3,2,2)
for iori = 1:4
    if ~isempty(int_all_ori{iori})
    hAll= cdfplot(int_all_ori{iori});
    hAll.Color = oriColMat(iori);
    hAll.LineStyle = '-';
    hold on
    end
end
xlabel('dF/F')
ylabel('cmlv frctn')
title('integral response - ori selective cells, all trials')
legend(oriStr)
print([fnout 'press_align_intOri_cdf' datasetStr '.pdf'], '-dpdf')

figure(allTrialsFig)
errAll_ori = cellfun(@(x) std(x,[],2)/sqrt(size(x,2)),resp_all_ori,'unif',false);
plt = 3:6;
for iori = 1:4
    subplot(3,2,plt(iori))
    shadedErrorBar(ttMs,mean(resp_all_ori{iori},2),errAll_ori{iori},'k')
    hold on
    xlim([-300 minTrialLengthFrames/(cycTime/cycTimeMs)])
    ylim([-0.01 0.05]);
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    title([oriStr(iori) ' slctv cells; n = ' num2str(size(resp_all_ori{iori},2))])
    xlabel('ms')
    ylabel('dF/F')
end

figure(allTrialsFig)
auc = zeros(nCells, noff, nDelta);
p_roc = zeros(nCells, noff, nDelta);
h_roc = zeros(nCells, noff, nDelta);
for ioff = 1:noff
    for itarg = 1:nDelta
        ind = intersect(find(tFramesOff(:,5) == offs(ioff)), find(targetDelta == deltas(itarg)));
        base = squeeze(bsxfun(@minus,mean(data_dfof(29:34,:,5,ind),1),mean(data_dfof(19:24,:,5,ind),1)));
        targ = squeeze(bsxfun(@minus,mean(data_dfof(29:34,:,6,ind),1),mean(data_dfof(19:24,:,6,ind),1)));
        for iC = 1:nCells
            auc(iC,ioff,itarg) = roc_gh(base(iC,:),targ(iC,:));
            [p_roc(iC,ioff,itarg) h_roc(iC,ioff,itarg)] = ranksum(base(iC,:),targ(iC,:));
        end
    end
end

figure;
auc_abs = auc-.5;
auc_norm = bsxfun(@rdivide, auc_abs, auc_abs(:,3,:)); 
for itarg = 1:nDelta
    subplot(2,1,itarg)
    for ioff = 1:noff
        errorbar(offs(ioff)*input.frameRateHz, mean(auc_norm(good_ind,ioff,itarg),1), std(auc_norm(good_ind,ioff,itarg),[],1)./sqrt(length(good_ind)), 'o')
        hold on
    end
    title([num2str(deltas(itarg)) ' deg target'])
    ylim([0 1])
end
suptitle('Base resp cells- normalized')

for itarg = 1:nDelta
    figure;
    for i = 1:length(good_ind)
        iC = good_ind(i);
        subplot(3,2,max_dir(iC,:))
        if sum(h_roc(iC,:,itarg),2)>0
            plot(auc_abs(iC,1,itarg),auc_abs(iC,3,itarg),'or')
        else
            plot(auc_abs(iC,1,itarg),auc_abs(iC,3,itarg),'ok')
        end
        hold on
    end
    for idir = 1:length(dirs)
        if length(find(max_dir(good_ind) == idir))>0
            subplot(3,2,idir)
            title(num2str(dirs(idir)))
            axis square
            plot(-0.4:0.01:0.4, -0.4:0.01:0.4, '-k')
            vline(0)
            hline(0)
            xlabel([num2str(offs(1)*input.frameRateHz) ' ms'])
            ylabel([num2str(offs(3)*input.frameRateHz) ' ms'])
        end
    end
    suptitle([num2str(deltas(itarg)) ' deg change'])
end

figure;
    for i = 1:length(good_ind)
        iC = good_ind(i);
        subplot(3,2,max_dir(iC,:))
        if sum(h_roc(iC,:,itarg),2)>0
            plot(auc_abs(iC,1,1),auc_abs(iC,1,2),'or')
        else
            plot(auc_abs(iC,1,1),auc_abs(iC,1,2),'ok')
        end
        hold on
    end
    for idir = 1:length(dirs)
        if length(find(max_dir(good_ind) == idir))>0
            subplot(3,2,idir)
            title(num2str(dirs(idir)))
            axis square
            plot(-0.4:0.01:0.4, -0.4:0.01:0.4, '-k')
            vline(0)
            hline(0)
            xlabel([num2str(deltas(1)) ' deg'])
            ylabel([num2str(deltas(2)) ' deg'])
        end
    end

    
    auc = zeros(nCells,noff);
for ioff = 1:noff
    base_list = [];
    basen1_list = [];
    for icyc = 3:nCyc-1
        ind = find(tFramesOff(:,icyc-1) == offs(ioff));
        base_list = cat(2,base_list,squeeze(mean(data_dfof_trial(29:34,:,icyc-1,ind),1)));
        basen1_list = cat(2,basen1_list,squeeze(mean(data_dfof_trial(29:34,:,icyc,ind),1)));
    end
    for iC = 1:nCells
        auc(iC,ioff) = roc_gh(base_list(iC,:), basen1_list(iC,:));
    end
end

auc = zeros(nCells,noff);
figure;
subplot(2,2,1)
for ioff = 1:noff
    base_list = [];
    basen1_list = [];
    for icyc = 3:nCyc
        ind = find(tFramesOff(:,icyc-1) == offs(ioff));
        base_list = cat(2,base_list,squeeze(mean(data_dfof_trial(29:34,:,icyc-1,ind),1)));
        basen1_list = cat(2,basen1_list,squeeze(mean(data_dfof_trial(29:34,:,icyc,ind),1)));
    end
    for iC = 1:nCells
        auc(iC,ioff) = roc_gh(base_list(iC,:), basen1_list(iC,:));
    end
end
errorbar(offs*input.frameRateHz,mean(auc(good_ind,:)-0.5,1), std(auc(good_ind,:)-0.5,[],1)./sqrt(length(good_ind)),'o')
title('All n-1 trials')
ylabel('AU-ROC')
xlabel('Interval (ms)')
ylim([-0.06 0.04])
auc = zeros(nCells,noff,noff);
for ioff1 = 1:noff
    subplot(2,2,ioff1+1)
    for ioff = 1:noff
        base_list = [];
        basen1_list = [];
        for icyc = 3:nCyc
            %ind = find(tFramesOff(:,icyc-1) == offs(ioff));
            ind = intersect(find(tFramesOff(:,icyc-1) == offs(ioff)),find(tFramesOff(:,icyc-2) == offs(ioff1)));
            if length(ind)>1
                base_list = cat(2,base_list,squeeze(mean(data_dfof_trial(29:34,:,icyc-1,ind),1)));
                basen1_list = cat(2,basen1_list,squeeze(mean(data_dfof_trial(29:34,:,icyc,ind),1)));
            else
                base_list = cat(2,base_list,squeeze(mean(data_dfof_trial(29:34,:,icyc-1,ind),1))');
                basen1_list = cat(2,basen1_list,squeeze(mean(data_dfof_trial(29:34,:,icyc,ind),1))');
            end
        end
        for iC = 1:nCells
            auc(iC,ioff,ioff1) = roc_gh(base_list(iC,:), basen1_list(iC,:));
        end
    end
    errorbar(offs*input.frameRateHz,mean(auc(ind0,:,ioff1)-0.5,1), std(auc(ind0,:,ioff1)-0.5,[],1)./sqrt(length(ind0)),'o')
    ylabel('AU-ROC')
    xlabel('Interval (ms)')
    title([num2str(offs(ioff1)*input.frameRateHz) ' n-1 trials'])
    ylim([-0.1 0.1])
end
suptitle('Baseline cycles 3-5')
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auroc_byInt.pdf']),'-dpdf')

auc = zeros(nCells,noff,nDelta);
for itarg = 1:nDelta
    for ioff1 = 1:noff
    for ioff = 1:noff
        ind = intersect(find(tFramesOff(:,4) == offs(ioff)), intersect(find(tFramesOff(:,5) == offs(ioff1)),find(targetDelta == deltas(itarg))));
        %ind = intersect(find(tFramesOff(:,5) == offs(ioff)),find(targetDelta == deltas(itarg)));
        base_list = squeeze(mean(data_dfof_trial(29:34,:,4,ind),1));
        targ_list = squeeze(mean(data_dfof_trial(29:34,:,5,ind),1));
        for i = 1:nCells
            auc(i,ioff,itarg) = roc_gh(base_list(i,:), targ_list(i,:));
        end
    end
    %good_ind_90 = intersect(find(max_dir == 4), good_ind);
    figure;
    subplot(2,2,1)
    errorbar(offs*input.frameRateHz,mean(auc(ind0,:,itarg)-0.5,1), (std(auc(ind0,:,itarg)-0.5,[],1))./sqrt(length(ind0)),'o')
    title(['ind0- ' num2str(length(ind0)) ' cells'])
    ylim([-0.2 0.2])
    subplot(2,2,2)
    errorbar(offs*input.frameRateHz,mean(auc(indint,:,itarg)-0.5,1), std(auc(indint,:,itarg)-0.5,[],1)./sqrt(length(indint)),'o')
    title(['indint- ' num2str(length(indint)) ' cells'])
    ylim([-0.2 0.2])
    subplot(2,2,3)
    errorbar(offs*input.frameRateHz,mean(auc(ind90,:,itarg)-0.5,1), std(auc(ind90,:,itarg)-0.5,[],1)./sqrt(length(ind90)),'o')
    title(['ind90- ' num2str(length(ind90)) ' cells'])
    ylim([-0.2 0.2])
    subplot(2,2,4)
    errorbar(offs*input.frameRateHz,mean(auc(good_resp_ind,:,itarg)-0.5,1), std(auc(good_resp_ind,:,itarg)-0.5,[],1)./sqrt(length(good_resp_ind)),'o')
    title(['All resp- ' num2str(length(good_resp_ind)) ' cells'])
    ylim([-0.2 0.2])
    suptitle([num2str(deltas(itarg)) ' deg change; n = ' num2str(offs(ioff1)*input.frameRateHz)])
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_auroc_byIntN1_' num2str(offs(ioff1)*input.frameRateHz) 'only_' num2str(deltas(itarg)) 'deg.pdf']),'-dpdf')
    end
end



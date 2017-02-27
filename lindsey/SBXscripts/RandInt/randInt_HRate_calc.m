
offT = nan(1,nTrials);
for itrial = 1:nTrials
    if ~FIx(itrial)
        temp = tFramesOff(itrial,find(~isnan(tFramesOff(itrial,:))));
        offT(1,itrial) = temp(tCyc(itrial));
    end
end

%hit rate calculation
hrate = zeros(noff,1);
ci95 = zeros(noff,2);
for ioff = 1:noff
    [hrate(ioff,:) ci95(ioff,:)] = binofit(sum(SIx(find(offT == offs(ioff))),2),(sum(SIx(find(offT == offs(ioff))),2) + sum(MIx(find(offT == offs(ioff))),2)));
end
figure; errorbar(offs.*input.frameRateHz, hrate', hrate-ci95(:,1), ci95(:,2)-hrate, '-ok')
ylim([0 1])

hrate_delt = zeros(nDelta,noff,1);
ci95_delt = zeros(nDelta,noff,2);
ind_n = zeros(2,nDelta,noff);
for ioff = 1:noff
    for itarg = 1:nDelta
        ind = intersect(find(targetDelta== deltas(itarg)),find(offT == offs(ioff)));
        [hrate_delt(itarg,ioff,:) ci95_delt(itarg,ioff,:)] = binofit(sum(SIx(ind),2),(sum(SIx(ind),2) + sum(MIx(ind),2)));
        ind_n(:,itarg,ioff) = [sum(SIx(ind),2);  sum(MIx(ind),2)];
    end
end
    
figure; 
subplot(1,2,1)
for itarg = 1:nDelta
    errorbar(offs.*input.frameRateHz, hrate_delt(itarg,:,:)', hrate_delt(itarg,:,:)-ci95_delt(itarg,:,1), ci95_delt(itarg,:,2)-hrate_delt(itarg,:,:), '-o')
    hold on
end
ylabel('Hit Rate')
xlabel('Interval (ms)')
ylim([0 1])
legend(num2str(deltas'), 'Location','Southeast')

subplot(1,2,2)
for ioff = 1:noff
    errorbar(deltas, hrate_delt(:,ioff,:), hrate_delt(:,ioff,:)-ci95_delt(:,ioff,1), ci95_delt(:,ioff,2)-hrate_delt(:,ioff,:), '-o')
    hold on
end
ylim([0 1])
ylabel('Hit Rate')
xlabel('Ori change (deg)')
suptitle([mouse ' ' date])
legend(num2str(offs.*input.frameRateHz), 'Location','Southeast')
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_hitRate_byInt.pdf']),'-dpdf', '-bestfit')

%FA rate calculation
cLeverUp = celleqel2mat_padded(input.cLeverUp);
CR = zeros(1,noff);
FA = zeros(1,noff);
FA_cyc = nan(1,nTrials);
FA_resp = cell(1,noff);
CR_resp = cell(1,noff);
for itrial = 1:nTrials
    if tCyc(itrial)>2
        tRelease = cLeverUp(itrial);
        for icyc = 3:tCyc(itrial)
            ind = find(offs == tFramesOff(itrial,icyc-1));
            if cStart(itrial)+sum(tFramesOff(itrial,1:icyc-1),2)+(input.nFramesOn{itrial}.*(icyc-1))+ceil(400./input.frameRateHz)>tRelease
                if cStart(itrial)+sum(tFramesOff(itrial,1:icyc-1),2)+(input.nFramesOn{itrial}.*(icyc-1))+ceil(100./input.frameRateHz)<tRelease
                    FA(1,ind) = FA(1,ind)+1;
                    FA_cyc(1,itrial) = icyc;
                    FA_resp{1,ind} = cat(3,FA_resp{1,ind}, data_dfof(:,:,icyc,itrial));
                end
            else
                CR(1,ind) = CR(1,ind)+1;
                CR_resp{1,ind} = cat(3,CR_resp{1,ind}, data_dfof(:,:,icyc,itrial));
            end
        end
    end
end

%FA rate by interval
[farate ci95] = binofit(FA,FA+CR);
figure; errorbar(offs*(input.frameRateHz), farate, farate'-ci95(:,1), ci95(:,2)+farate', 'o')
ylim([0 1])
ylabel('FA Rate')
xlabel('Interval (ms)')
suptitle([mouse ' ' date])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FARate_byInt.pdf']),'-dpdf', '-bestfit')

%plot TCs of FAs vs CRs
figure;
for ioff = 1:noff
    subplot(2,1,1)
    temp = mean(FA_resp{1,ioff},3);
    plot(tt,mean(bsxfun(@minus,temp(:,good_ind),mean(temp(base_win,good_ind),1)),2))
    ylim([-.02 0.05])
    hold on
    subplot(2,1,2)
    temp = mean(CR_resp{1,ioff},3);
    plot(tt,mean(bsxfun(@minus,temp(:,good_ind),mean(temp(base_win,good_ind),1)),2))    
    hold on
    ylim([-.02 0.05])
end
subplot(2,1,1)
vline(([resp_win(1) resp_win(end)]-20)*frameRateHz)
xlabel('Time from stimulus')
ylabel('dF/F')
title(['False Alarms- ' num2str(FA)])
legend(num2str(offs.*input.frameRateHz), 'Location','Northwest')
subplot(2,1,2)
vline(([resp_win(1) resp_win(end)]-20)*frameRateHz)
xlabel('Time from stimulus')
ylabel('dF/F')
title(['Correct Rejects- ' num2str(CR)])
legend(num2str(offs.*input.frameRateHz), 'Location','Northwest')
suptitle([mouse ' ' date])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TC_byFA_byInt.pdf']),'-dpdf', '-bestfit')

Miss = zeros(nDelta,noff);
Hit = zeros(nDelta,noff);
RT = cell(nDelta,noff);
Miss_resp = cell(nDelta,noff);
Hit_resp = cell(nDelta,noff);
for itrial = 1:nTrials
    ind = find(offs == tFramesOff(itrial,tCyc(itrial)));
    itarg = find(deltas == targetDelta(itrial));
    if SIx(itrial)
        if input.reactTimesMs{itrial}>300
            RT{itarg,ind} = [RT{itarg,ind} input.reactTimesMs{itrial}];
            Hit(itarg,ind) = Hit(itarg,ind)+1;
            Hit_resp{itarg,ind} = cat(3,Hit_resp{itarg,ind}, data_dfof(:,:,nCyc(itrial),itrial));
        end
    elseif MIx(itrial)
        Miss(itarg,ind) = Miss(itarg,ind)+1;
        Miss_resp{itarg,ind} = cat(3,Miss_resp{itarg,ind}, data_dfof(:,:,nCyc(itrial),itrial));
    end
end
figure;
for itarg = 1:nDelta
    for ioff = 1:noff
        if Hit(itarg,ioff)>0
            subplot(2,2,1+((itarg-1)*nDelta))
            temp = mean(Hit_resp{itarg,ioff},3);
            plot(tt,mean(bsxfun(@minus,temp(:,good_resp_ind),mean(temp(base_win,good_resp_ind),1)),2))
            ylim([-.02 0.075])
            hold on
        end
        if Miss(itarg,ioff)>0
            subplot(2,2,2+((itarg-1)*nDelta))
            temp = mean(Miss_resp{itarg,ioff},3);
            plot(tt,mean(bsxfun(@minus,temp(:,good_resp_ind),mean(temp(base_win,good_resp_ind),1)),2))    
            hold on
            ylim([-.02 0.075])
        end
    end
    subplot(2,2,1+((itarg-1)*nDelta))
    title(['Hits- ' num2str(deltas(itarg)) ' deg- ' num2str(Hit(itarg,:))])
    vline(([resp_win(1) resp_win(end)]-20)*frameRateHz)
    xlabel('Time from target')
    ylabel('dF/F')
    legend(num2str(offs.*input.frameRateHz), 'Location','Northwest')
    subplot(2,2,2+((itarg-1)*nDelta))
    title(['Misses- ' num2str(deltas(itarg)) ' deg- ' num2str(Miss(itarg,:))])
    vline(([resp_win(1) resp_win(end)]-20)*frameRateHz)
    xlabel('Time from target')
    ylabel('dF/F')
    legend(num2str(offs.*input.frameRateHz), 'Location','Northwest')
end
suptitle([mouse ' ' date])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TC_byHit_byInt.pdf']),'-dpdf', '-bestfit')


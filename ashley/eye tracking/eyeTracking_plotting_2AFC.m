tLeftTrial = cell2mat(input.tLeftTrial);
tLeftResponse  = cell2mat(input.tLeftResponse);
tRightResponse = cell2mat(input.tRightResponse);
centroid_mat_start = permute(centroid_mat_start,[1,3,2]);
centroid_mat_decide = permute(centroid_mat_decide,[1,3,2]);

rad_mat_norm = rad_mat_start./nanmean(nanmean(rad_mat_start(preevent_frames/2:preevent_frames,:),1),2);
centroid_mat_norm = bsxfun(@minus,centroid_mat_start,nanmean(nanmean(centroid_mat_start(preevent_frames/2:preevent_frames,:,:),1),2));
centroid_mat_norm = permute(centroid_mat_norm,[1,3,2]);

indLcorr = intersect(find(tLeftResponse),find(tLeftTrial));
indLincorr = intersect(find(tLeftResponse),find(~tLeftTrial));
indRcorr = intersect(find(tRightResponse),find(~tLeftTrial));
indRincorr = intersect(find(tRightResponse),find(tLeftTrial));

figure;
subplot(3,2,1)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLincorr),2)', nanstd(rad_mat_norm(:,indLincorr),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLcorr),2)', nanstd(rad_mat_norm(:,indLcorr),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,2)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRincorr),2)', nanstd(rad_mat_norm(:,indRincorr),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRcorr),2)', nanstd(rad_mat_norm(:,indRcorr),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,3)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLincorr)),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLcorr)),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,4)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRincorr)),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRcorr)),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,5)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLincorr)),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLcorr)),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,6)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRincorr)),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRcorr)),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])
suptitle([mouse ' ' date '- trial start align'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_byResp.pdf']),'-dpdf', '-bestfit')

indLcorr = intersect(find(tLeftResponse),find(tLeftTrial));
indLincorr = intersect(find(tLeftResponse),find(~tLeftTrial));
indRcorr = intersect(find(tRightResponse),find(~tLeftTrial));
indRincorr = intersect(find(tRightResponse),find(tLeftTrial));

figure;
subplot(3,2,1)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRcorr),2)', nanstd(rad_mat_norm(:,indRcorr),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLcorr),2)', nanstd(rad_mat_norm(:,indLcorr),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,2)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLincorr),2)', nanstd(rad_mat_norm(:,indLincorr),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRincorr),2)', nanstd(rad_mat_norm(:,indRincorr),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,3)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRcorr)),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLcorr)),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,4)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLincorr)),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRincorr)),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,5)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRcorr)),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLcorr)),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,6)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLincorr)),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRincorr)),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])
suptitle([mouse ' ' date '- trial start align'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_bySide.pdf']),'-dpdf', '-bestfit')

probList = cell2mat(input.ProbList);
tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
col_mat = strvcat('b', 'k', 'g');
indL_n = [];
indR_n = [];
figure;
for iprob = 1:length(probList)
    prob = probList(iprob);
    ind = find(tProbLeft == prob);
    ind = 1+((iprob-1)*80):80+((iprob-1)*80);
    indL = intersect(ind, intersect(find(tLeftResponse),find(tLeftTrial)));
    indR = intersect(ind, intersect(find(tRightResponse),find(~tLeftTrial)));
    indL_n = [indL_n length(indL)];
    indR_n = [indR_n length(indR)];
    subplot(3,2,1)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(rad_mat_norm(:,indL),2)', nanstd(rad_mat_norm(:,indL),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,2)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(rad_mat_norm(:,indR),2)', nanstd(rad_mat_norm(:,indR),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
    subplot(3,2,3)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,1,indL)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indL)),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,4)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,1,indR)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indR)),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
    subplot(3,2,5)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,2,indL)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indL)),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,6)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,2,indR)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indR)),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
end
subplot(3,2,1)
title(['Left trials- ' num2str(indL_n)])
ylim([.95 1.25])
ylabel('Radius')
subplot(3,2,2)
title(['Right trials- ' num2str(indR_n)])
ylim([.95 1.25])    
ylabel('Radius')
subplot(3,2,3)
ylim([-.1 .2])
ylabel('Horizontal Position')
subplot(3,2,4)
ylim([-.1 .2])
ylabel('Horizontal Position')
subplot(3,2,5)
ylim([-.1 .2])
ylabel('Vertical Position')
subplot(3,2,6)
ylim([-.1 .2])
ylabel('Vertical Position')
suptitle([mouse ' ' date '- trial start align; Block order: black, blue, cyan, green'])
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_byBlock.pdf']),'-dpdf', '-bestfit')

edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m
edit successTrials.m

dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_ind = cycV_ind{1};
AV_ind = cycAV_ind{1};
Success_ind = cycSuccess_ind{1};

Success_ind = find(strcmp('success',trialOutcome));
Miss_ind = find(strcmp('ignore',trialOutcome));
Early_ind = find(strcmp('failure',trialOutcome));
%%

V_Success_ind = find(ismember(Success_ind,intersect(Success_ind,V_ind)));
AV_Success_ind = find(ismember(Success_ind,intersect(Success_ind,AV_ind)));

V_successAvg = mean(mean(dataTrialStart(:,:,V_Success_ind),3),2);
AV_successAvg = mean(mean(dataTrialStart(:,:,AV_Success_ind),3),2);

errbar_V = std(mean(dataTrialStart(:,:,V_Success_ind),2),[],3)/sqrt(size(dataTrialStart(:,:,V_Success_ind),3));
errbar_AV = std(mean(dataTrialStart(:,:,AV_Success_ind),2),[],3)/sqrt(size(dataTrialStart(:,:,AV_Success_ind),3));

figure;
% plot(V_successAvg,'g');
shadedErrorBar([],V_successAvg,errbar_V,'g', 1);
hold on
% plot(AV_successAvg,'m');
shadedErrorBar([],AV_successAvg,errbar_AV,'m', 1);
hold on
vline(30,'k')
hold on
vline(30+cycTime,'k:');
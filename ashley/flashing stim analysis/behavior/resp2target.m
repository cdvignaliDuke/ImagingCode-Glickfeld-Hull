% edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m
edit FlashingStim_dataSortedByCycle_combineDatasets.m


Success_ind = find(strcmp('success',trialOutcome));
Miss_ind = find(strcmp('ignore',trialOutcome));
Early_ind = find(strcmp('failure',trialOutcome));

Dirs = unique(DirectionDeg);

%%
L = 60;
data_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
dF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind)); 
dFoverF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
for i = 1:length(Success_ind)
    trial = Success_ind(i);
    data_aroundTarget(:,:,i) = dataTC(cTargetOn(trial)-(L/2):cTargetOn(trial)+((L/2)-1),:);
    dF_aroundTarget(:,:,i) = bsxfun(@minus, data_aroundTarget(:,:,i), mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
    dFoverF_aroundTarget(:,:,i) = bsxfun(@rdivide, dF_aroundTarget(:,:,i),mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
end

V_Success_ind = find(ismember(Success_ind,intersect(Success_ind,V_ind)));
AV_Success_ind = find(ismember(Success_ind,intersect(Success_ind,AV_ind)));
% A_Success_ind = find(ismember(Success_ind,intersect(Success_ind,A_ind)));

V_successAvg = mean(mean(dFoverF_aroundTarget(:,:,V_Success_ind),3),2);
AV_successAvg = mean(mean(dFoverF_aroundTarget(:,:,AV_Success_ind),3),2);
% A_successAvg = mean(mean(dFoverF_aroundTarget(:,:,A_Success_ind),3),2);

errbar_V = std(mean(dFoverF_aroundTarget(:,:,V_Success_ind),2),[],3)/sqrt(size(dFoverF_aroundTarget(:,:,V_Success_ind),3));
errbar_AV = std(mean(dFoverF_aroundTarget(:,:,AV_Success_ind),2),[],3)/sqrt(size(dFoverF_aroundTarget(:,:,AV_Success_ind),3));
% errbar_A = std(mean(dFoverF_aroundTarget(:,:,A_Success_ind),2),[],3)/sqrt(size(dFoverF_aroundTarget(:,:,A_Success_ind),3));

figure;
% plot(V_successAvg,'g');
% plot(AV_successAvg,'m');
% errorbar(A_successAvg,errbar_A,'r');
% hold on
errorbar(AV_successAvg,errbar_AV,'m');
hold on
errorbar(V_successAvg,errbar_V,'g');
hold on
vline((L/2),'c')
hold on
for i = 1:2
    vline((L/2)-(cycTime*i),'k:');
end
hold on
% vline((L/2)+tooFastTime,'k')
% vline((L/2) + maxReactTime,'k')
vline((L/2) + meanSuccessReactTime, 'b--')
hold on
ylabel('dF/F')
xlabel('frames')
title('average target response, all cells')



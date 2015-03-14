edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m

Success_ind = find(strcmp('success',trialOutcome));
Miss_ind = find(strcmp('ignore',trialOutcome));
Early_ind = find(strcmp('failure',trialOutcome));



figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = intersect(find(tCyclesOn == cycles(icyc)),Success_ind);
    V_cycInd = find(ismember(trials,V_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,:,V_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' auditory trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

for icyc = 1:length(cycles)
    dataDFoverF_cmlvNoTarget = [];
    V_indAll = [];
    AV_indAll = [];
    Success_indAll = [];
    running_ind = 0;
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        dataDFoverF = cycDataDFoverF{i};
        dataDFoverF_NoTarget = zeros(L,size(dataDFoverF,2),size(dataDFoverF,3));
        dataDFoverF_NoTarget = dataDFoverF(1:L,:,:);
        dataDFoverF_cmlvNoTarget = cat(3,dataDFoverF_cmlvNoTarget,dataDFoverF_NoTarget);
        trials = find(tCyclesOn == cycles(i));
        V_cycInd = find(ismember(trials,V_ind));
        AV_cycInd = find(ismember(trials,AV_ind));
        Success_cycInd = find(ismember(trials,Success_ind));
        V_indAll = cat(1,V_indAll, V_cycInd+running_ind);
        AV_indAll = cat(1,AV_indAll,AV_cycInd+running_ind);
        Success_indAll = cat(1,Success_indAll,Success_cycInd+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycAV_ind{icyc} = AV_indAll;
    cycSuccess_ind{icyc} = Success_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    Success_cycInd = cycSuccess_ind{icyc};
    Vind = intersect(V_cycInd,Success_cycInd);
    AVind = intersect(AV_cycInd,Success_cycInd);
    V_avg = mean(mean(dataDFoverF(:,:,Vind),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AVind),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(AV_avg,'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' auditory trials'])
    hold on
end

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

V_successAvg = mean(mean(dFoverF_aroundTarget(:,:,V_Success_ind),3),2);
AV_successAvg = mean(mean(dFoverF_aroundTarget(:,:,AV_Success_ind),3),2);

figure;
plot(V_successAvg,'g');
hold on
plot(AV_successAvg,'m');
hold on
vline((L/2),'c')
hold on
for i = 1:2
    vline((L/2)-(cycTime*i),'k:');
end
hold on
vline((L/2)+tooFastTime,'k')
vline((L/2) + maxReactTime,'k')
vline((L/2) + meanSuccessReactTime, 'b')
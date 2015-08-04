edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m

Success_ind = find(strcmp('success',trialOutcome));
Miss_ind = find(strcmp('ignore',trialOutcome));
Early_ind = find(strcmp('failure',trialOutcome));



figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = intersect(find(tCyclesOn == cycles(icyc)),Early_ind);
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
    Early_indAll = [];
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
        Early_cycInd = find(ismember(trials,Early_ind));
        V_indAll = cat(1,V_indAll, V_cycInd+running_ind);
        AV_indAll = cat(1,AV_indAll,AV_cycInd+running_ind);
        Early_indAll = cat(1,Early_indAll,Early_ind+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycAV_ind{icyc} = AV_indAll;
    cycEarly_ind{icyc} = Early_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    Early_cycInd = cycEarly_ind{icyc};
    Vind = intersect(V_cycInd,Early_cycInd);
    AVind = intersect(AV_cycInd,Early_cycInd);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(Vind)) ' visual trials; ' num2str(length(AVind)) ' auditory trials'])
    hold on
end
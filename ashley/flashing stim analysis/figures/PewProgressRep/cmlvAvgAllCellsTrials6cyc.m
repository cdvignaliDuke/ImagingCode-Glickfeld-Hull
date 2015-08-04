cycleLength = 6;
cyc = find(cycles == cycleLength);
dataCellsTrials = cycDataDFoverF_cmlvNoTarget{cyc};
V_cycInd = cycV_ind{cyc};
A_cycInd = cycA_ind{cyc};
AV_cycInd = cycAV_ind{cyc};
V_avg = mean(mean(dataCellsTrials(:,:,V_cycInd),3),2);
A_avg = mean(mean(dataCellsTrials(:,:,A_cycInd),3),2);
AV_avg = mean(mean(dataCellsTrials(:,:,AV_cycInd),3),2);
errbar_V = (std(mean(dataCellsTrials(:,:,V_cycInd),2),[],3))/(sqrt(size(dataCellsTrials(:,:,V_cycInd),3)));
errbar_A = (std(mean(dataCellsTrials(:,:,A_cycInd),2),[],3))/(sqrt(size(dataCellsTrials(:,:,A_cycInd),3)));
errbar_AV = (std(mean(dataCellsTrials(:,:,AV_cycInd),2),[],3))/(sqrt(size(dataCellsTrials(:,:,AV_cycInd),3)));
legendinfo = {'vis only','aud only','vis+aud'};

figure;
shadedErrorBar([],V_avg(20:end,:),errbar_V(20:end,:),'k',1)
hold on
shadedErrorBar([],A_avg(20:end,:),errbar_A(20:end,:),'b',1)
hold on
shadedErrorBar([],AV_avg(20:end,:),errbar_AV(20:end,:),'g',1)
hold on
vline(10,'k')
hold on
for i = 1:cycles(icyc)-1
    L = (i*cycTime)+11;
    vline(L,'k:');
    hold on
end
vline((cycles(icyc)*cycTime+11),'c');
hold on
title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(A_cycInd)) ' auditory trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
hold on
xlim([0 length(V_avg(20:end,:))+5])
ylim([-0.05 0.05])
if icyc == 1
    legend(legendinfo,'Location','NorthWest')
end
group = intersect(baselineStimRespIndex_V,cellsPrefRespZero);
for i = 1:length(group)
cell = group(i);
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,cell,V_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cell,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,cell,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cell,V_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,cell,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cell,AV_cycInd),3)));
    subplot(3,4,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
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
    title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])

end
end

%%
%success only trials

group = intersect(baselineStimRespIndex_V,cellsPrefRespZero);
for i = 1:length(group)
cell = group(i);
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    V_avg = mean(mean(dataDFoverF(:,cell,V_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cell,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,cell,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cell,V_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,cell,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cell,AV_cycInd),3)));
    subplot(3,4,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
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
    title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])

end
end
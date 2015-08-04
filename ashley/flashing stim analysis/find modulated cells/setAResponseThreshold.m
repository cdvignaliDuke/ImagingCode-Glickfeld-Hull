figure;
icyc = 10;
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,:,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,V_cycInd),3)));
%     errbar_A = (std(mean(dataDFoverF(:,:,A_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,A_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,:,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,AV_cycInd),3)));
%     subplot(3,4,icyc);
%     plot(V_avg(20:end,:),'g');
%     hold on
%     plot(A_avg(20:end,:),'r');
%     hold on
%     plot(AV_avg(20:end,:),'m');
%     hold on
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
%     errorbar(A_avg(20:end,:),errbar_A(20:end,:),'b')
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
    if icyc == 1
        title([num2str(size(dataDFoverF,2)) ' cells'])
    else
%     title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(A_cycInd)) ' auditory trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
%     if icyc == 1
%         legend(legendinfo,'Location','SouthEast')
%     end

%%
%plot each trial
figure;
datasmooth1 = tsmovavg(squeeze(mean(dataDFoverF(:,baselineStimRespIndex_V,V_cycInd),3)),'s',3,1);
dataavgtrials = squeeze(mean(dataDFoverF(:,:,V_cycInd),3));
plot(dataavgtrials)
vline(30,'k') 
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+31;
        vline(L,'k:');
        hold on
    end
    
drivencells = find(any(dataavgtrials >0.05,1));

figure;
plot(dataavgtrials(:,drivencells(2)))
vline(30,'k') 
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+31;
        vline(L,'k:');
        hold on
    end
%%   
for i = 1:length(drivencells)
figure;
% 
%     subplot(5,5,i)
plot(tsmovavg(squeeze(dataDFoverF(:,(drivencells(i)),V_cycInd)),'s',3,1),'g')
alpha(0.25)
hold on
plot(tsmovavg(squeeze(dataDFoverF(:,(drivencells(i)),AV_cycInd)),'s',3,1),'k')
% hold on
vline(30,'k') 
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+31;
        vline(L,'k:');
        hold on
    end
end
%%
figure;
hold on
plot(mean(tsmovavg(squeeze(dataDFoverF(:,(drivencells(38)),V_cycInd)),'s',3,1),2),'r')
hold on
plot(mean(tsmovavg(squeeze(dataDFoverF(:,(drivencells(38)),AV_cycInd)),'s',3,1),2),'r')

    
    
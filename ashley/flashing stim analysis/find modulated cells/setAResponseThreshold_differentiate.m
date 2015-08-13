figure;
icyc = 7;
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
 
for i = 1:length(baselineStimRespIndex_V)
figure;
% subplot(2,1,i)
    
plot(tsmovavg(squeeze(dataDFoverF(:,(baselineStimRespIndex_V(i)),V_cycInd)),'s',3,1),'g','LineWidth',2)
hold on
plot(tsmovavg(squeeze(dataDFoverF(:,(baselineStimRespIndex_V(i)),AV_cycInd)),'s',3,1),'k')
% hold on
vline(30,'k') 
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+31;
        vline(L,'k:');
        hold on
    end
    xlim([1 size(dataDFoverF,1)])
end

%%
datadiff = zeros(size(dataDFoverF));
for icell = 1:size(dataDFoverF,2)
for itrial = 1:size(dataDFoverF,3)
    datadiff(2:end,icell,itrial) = diff(squeeze(dataDFoverF(:,icell,itrial)));
end
end
datadiff(datadiff < 0.05) = 0;

%%

cell1 = 3;

for i = 1:size(dataDFoverF,3)
if ismember(i,[1:15:nCells]) == 1
    figure;
    plotN = 1;
end
    subplot(5,3,plotN)
    
plot(tsmovavg(squeeze(dataDFoverF(:,cell1,i)),'s',3,1),'k')
hold on
plot(squeeze(dataDFoverF(:,cell1,i)),'b')
hold on
plot(squeeze(datadiff(:,cell1,i)),'r')
hline(0,'k:')
vline(30,'k') 
xlim([1 size(dataDFoverF,1)])
title(['trial ' num2str(i) '; cell ' num2str(cell1)])
plotN = plotN+1;
end






%% differentiate response for each trial
data = cycData{icyc};
dataDF = cycDataDF{icyc};
datadiff = diff(data);
datadiff_norm = diff(dataDFoverF);

data_spike = datadiff_norm >= 0.1;

%%
data_spike = datadiff_norm >= 0.1;
cell1 = 90;
data_spikes_cell1 = squeeze(data_spike(:,cell1,i));

figure;
for j = 1:length(drivencells)
    
for i = 1:size(data_spikes_cell1,2);
    spike_fr = find(squeeze(data_spike(:,j,i)) == 1);
    subplot(4,5,j)
    if isempty(spike_fr) == 0
        if isempty(find(V_cycInd == i)) == 0
            
            plot(spike_fr,i,'go')
        else

            plot(spike_fr,i,'ko')
        end
        hold on
        vline(30,'k') 
        for i = 1:cycles(icyc)-1
            L = (i*cycTime)+31;
            vline(L,'k:');
            hold on
        end
    end
xlim([1 size(dataDFoverF,1)])    
ylim([1 size(dataDFoverF,3)])
end
end


%% try downsampling

squeeze(mean(reshape(timecourse, [5 siz(1)/5 siz(2)]),1))

%%
    
figure;
for i = 1:size(data,3)
    subplot(5,5,i)
plot(tsmovavg(squeeze(dataDFoverF(:,drivencells(4),i)),'s',3,1),'g');
hold on
plot(datadiff_norm(:,drivencells(4),i),'r')
hold on
vline(30,'k') 
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+31;
        vline(L,'k:');
        hold on
    end
    xlim([1 size(dataDFoverF,1)])
end



    
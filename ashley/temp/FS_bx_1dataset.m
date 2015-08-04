SubNum = '613';
date = '150511';
time = '1453';
ImgFolder = '003';
mouse = 'AW13';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% load dataTC
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('dataTC.mat');
load('Timecourses.mat');
% dataTC = dataTimecourse.dataTC - dataTimecourse.npilTC ;
% dataTC(dataTC<0) = 0;
dataTC = dataTimecourse.dataTC;


%variables from mworks
cLeverDown = cell2mat_padded(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
tCyclesOn = cell2mat_padded(input.tCyclesOn);
block2 = cell2mat_padded(input.tBlock2TrialNumber);
cycles = unique(tCyclesOn);
cycTime = input.nFramesOn + input.nFramesOff;
% frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
% RateFRperMS = frameRateS/1000;
% cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;
trialOutcome = input.trialOutcomeCell;

% % special case where last trial take place within last 2 frames collected
% % by mworks, but not scanbox
% cLeverDown = cLeverDown(1:end-2,:);
% cTargetOn = cTargetOn(1:end-2,:);
% tCyclesOn = tCyclesOn(1:end-2,:);
% block2 = block2(1:end-2,:);
% nTrials = nTrials-2;

%dFoverF by cycle length
for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    Data = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDFoverFAll = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    for itrial = 1:length(ind)
        Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverFAll(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
    cycData{icyc} = Data;
    cycDataDF{icyc} = DataDF;
    cycDataDFoverF{icyc} = DataDFoverFAll;
end

%% find successful trials
for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    suc_ind = strcmp(trialOutcome(ind),'success');
    data = cycDataDFoverF{icyc};
    cycDataDFoverF_success{icyc} = data(:,:,suc_ind);
    V_ind_success{icyc} = find(block2(ind(suc_ind)) == 0);
    A_ind_success{icyc} = find(block2(ind(suc_ind)) == 1);
end

%%
% plot average of all cells for like-cycle lengths
figure;
for icyc = 1:length(cycles)
    DataDFoverFAll = cycDataDFoverF_success{icyc};
    if isempty(DataDFoverFAll)
        subplot(3,3,icyc);
        hline(0,'k')
    else
    V_ind = V_ind_success{icyc};
    A_ind = A_ind_success{icyc};
    V_avg = mean(mean(DataDFoverFAll(:,:,V_ind),3),2);
    A_avg = mean(mean(DataDFoverFAll(:,:,A_ind),3),2);
    V_cellsAvg = squeeze(mean(DataDFoverFAll(:,:,V_ind),2));
    A_cellsAvg = squeeze(mean(DataDFoverFAll(:,:,A_ind),2));
    subplot(3,3,icyc);
    plot(V_avg,'g','LineWidth',3);
    hold on
%     plot(V_cellsAvg, 'g');
    hold on
    plot(A_avg,'m','LineWidth',3);
    hold on
%     plot(A_cellsAvg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials']);
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
end

% plot cumulative average for each trial length (# of cycles), average
% up to, but not including the target

L = zeros(size(cycles));
for icyc = 1:length(cycles)
    DataDFoverFAll = [];
    V_indAll = [];
    A_indAll = [];
    running_ind = 0;
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        dataDFoverF = cycDataDFoverF_success{i};
        dataDFoverF_NoTarget = zeros(L,size(dataDFoverF,2),size(dataDFoverF,3));
        dataDFoverF_NoTarget = dataDFoverF(1:L,:,:);
        DataDFoverFAll = cat(3,DataDFoverFAll,dataDFoverF_NoTarget);
        V_ind = V_ind_success{i};
        A_ind = A_ind_success{i};
        trials = length(V_ind) + length(A_ind);
        V_indAll = cat(1,V_indAll, V_ind+running_ind);
        A_indAll = cat(1,A_indAll, A_ind+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget_success{icyc} = DataDFoverFAll; 
    cycV_ind_success{icyc} = V_indAll;
    cycA_ind_success{icyc} = A_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget_success{icyc};
    if isempty(dataDFoverF)
        subplot(3,3,icyc);
        hline(0,'k');
    else
    V_ind = cycV_ind_success{icyc};
    A_ind = cycA_ind_success{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_ind),3),2);
    A_avg = mean(mean(dataDFoverF(:,:,A_ind),3),2);
    V_cellsAvg = squeeze(mean(dataDFoverF(:,:,V_ind),2));
    A_cellsAvg = squeeze(mean(dataDFoverF(:,:,A_ind),2));
    subplot(3,3,icyc);
    hold on
%     plot(V_cellsAvg, 'g')
    hold on
    plot(V_avg,'g','LineWidth',3);
    hold on
%     plot(A_cellsAvg,'r')
    hold on
    plot(A_avg,'m','LineWidth', 3);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
    hold on
    end
end

%% use error bars
figure;
V_ind = find(block2 == 0);
AV_ind = find(block2 == 1);
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,:,V_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,:,AV_cycInd),3),2);
    errbar_V = (std(mean(DataDFoverF(:,:,V_cycInd),2),[],3))/(sqrt(size(DataDFoverF(:,:,V_cycInd),3)));
    errbar_AV = (std(mean(DataDFoverF(:,:,AV_cycInd),2),[],3))/(sqrt(size(DataDFoverF(:,:,AV_cycInd),3)));
    subplot(3,3,icyc);
    errorbar(V_avg,errbar_V,'g');
    hold on
    errorbar(AV_avg,errbar_AV,'m');
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

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    errbar_V = std(mean(dataDFoverF(:,:,V_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,:,V_cycInd),3));
%     errbar_A = std(mean(dataDFoverF(:,cellsPrefZero,A_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefZero,A_cycInd),3));
    errbar_AV = std(mean(dataDFoverF(:,:,AV_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,:,AV_cycInd),3));
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,cellsPrefZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,:,AV_cycInd),3),2);
    subplot(3,3,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g');
    hold on
%     errorbar(A_avg(20:end,:),errbar_A(20:end,:),'r');
%     hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'m');
    hold on
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+10;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+10),'c');
    hold on
%     title(['n = ' num2str(length(cellsPrefZero)) '; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
%     title(['n = ' num2str(size(dataDFoverF,2) '; ' num2str(length(V_cycInd)) ' vis trials; '  num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end

% success only*****
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget_success{icyc};
    if isempty(dataDFoverF)
        subplot(3,3,icyc);
        hline(0,'k');
    else
    V_ind = cycV_ind_success{icyc};
    A_ind = cycA_ind_success{icyc};
    errbar_V = (std(mean(dataDFoverF(:,:,V_ind),2),[],3))/(sqrt(size(dataDFoverF(:,:,V_ind),3)));
    errbar_A = (std(mean(dataDFoverF(:,:,A_ind),2),[],3))/(sqrt(size(dataDFoverF(:,:,A_ind),3)));
    V_avg = mean(mean(dataDFoverF(:,:,V_ind),3),2);
    A_avg = mean(mean(dataDFoverF(:,:,A_ind),3),2);
    subplot(3,3,icyc);
    errorbar(V_avg,errbar_V,'g')
    hold on
    errorbar(A_avg,errbar_A,'m')
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
    hold on
    end
end


%%
cycTrialLength = 7;
cyc = find(cycles == cycTrialLength);
dataCellsTrials = cycDataDFoverF{cyc};
V_ind = V_ind_success{cyc};
A_ind = A_ind_success{cyc};


cellsperpage = 16;
pages = ceil(size(dataCellsTrials,2)/cellsperpage);


cell = 1;
for ifig = 1:pages
figure;
for iplot = 1:cellsperpage
    subplot(4,4,iplot)
%     plot(squeeze(dataCellsTrials(:,cell,V_ind)),'g');
%     hold on
    plot(mean(dataCellsTrials(:,cell,V_ind),3),'g','LineWidth',3)
    hold on
%     plot(squeeze(dataCellsTrials(:,cell,A_ind)),'m')
%     hold on
    plot(mean(dataCellsTrials(:,cell,A_ind),3),'m','LineWidth',3)
    hold on
    vline(30,'k');
    hold on
    for i = 1:cyc-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cyc*cycTime+30),'c');
    xlim([0 size(dataCellsTrials,1)])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_ind)) ' vis, ' num2str(length(A_ind)) 'aud'])
    end
%     if ismember(cell,cellsPrefZero) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
    cell = cell+1;
end    
end
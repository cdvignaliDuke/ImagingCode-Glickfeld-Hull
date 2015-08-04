SubNum = '618';
date = '150526';
time = '1520';
ImgFolder = '004';
mouse = 'AW18';

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
% dataTC = dataTimecourse.dataTC;
dataTC = dataTimecourse.dataTCsub;


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
    DataDFoverF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    for itrial = 1:length(ind)
        Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
    cycData{icyc} = Data;
    cycDataDF{icyc} = DataDF;
    cycDataDFoverF{icyc} = DataDFoverF;
end

% save('cycDataDFoverF.mat','cycDataDFoverF');

% plot average of all cells for like-cycle lengths
figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_ind = find(block2(trials) == 0);
    A_ind = find(block2(trials) == 1);
    V_avg = mean(mean(DataDFoverF(:,:,V_ind),3),2);
    A_avg = mean(mean(DataDFoverF(:,:,A_ind),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
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

% plot cumulative average for each trial length (# of cycles), average
% up to, but not including the target

L = zeros(size(cycles));
for icyc = 1:length(cycles)
    dataDFoverF_cmlvNoTarget = [];
    V_indAll = [];
    A_indAll = [];
    running_ind = 0;
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        dataDFoverF = cycDataDFoverF{i};
        dataDFoverF_NoTarget = zeros(L,size(dataDFoverF,2),size(dataDFoverF,3));
        dataDFoverF_NoTarget = dataDFoverF(1:L,:,:);
        dataDFoverF_cmlvNoTarget = cat(3,dataDFoverF_cmlvNoTarget,dataDFoverF_NoTarget);
        trials = find(tCyclesOn == cycles(i));
        V_ind = find(block2(trials) == 0);
        A_ind = find(block2(trials) == 1);
        V_indAll = cat(1,V_indAll, V_ind+running_ind);
        A_indAll = cat(1,A_indAll, A_ind+running_ind);
        running_ind = length(trials)+running_ind;
    end        
    cycDataDFoverF_cmlvNoTarget{icyc} = dataDFoverF_cmlvNoTarget; 
    cycV_ind{icyc} = V_indAll;
    cycA_ind{icyc} = A_indAll;
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_ind = cycV_ind{icyc};
    A_ind = cycA_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_ind),3),2);
    A_avg = mean(mean(dataDFoverF(:,:,A_ind),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
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

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,:,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,:,A_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,:,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,V_cycInd),3)));
    errbar_A = (std(mean(dataDFoverF(:,:,A_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,A_cycInd),3)));
    subplot(3,3,icyc);
%     plot(V_avg(20:end,:),'g');
%     hold on
%     plot(A_avg(20:end,:),'r');
%     hold on
    shadedErrorBar([],V_avg,errbar_V,'g',1)
    hold on
    shadedErrorBar([],A_avg,errbar_A,'k',1)
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
        title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(A_cycInd)) ' aud trials'])
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end

% save('cycDataDFoverF_cmlvNoTarget.mat', 'cycDataDFoverF_cmlvNoTarget', 'cycA_ind', 'cycV_ind');

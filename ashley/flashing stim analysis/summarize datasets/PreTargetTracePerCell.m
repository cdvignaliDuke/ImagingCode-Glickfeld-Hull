clear all
close all
% fnouttemp = ['Z:\Analysis\temp figs\150924retreatposter']; %
awFSAVdatasets
% for iexp = 1:size(expt,2)
iexp = 8;
    
%%
ialign = 1;
run('divideupdatabyalignment.m')

%%
%%
analysisName = '';
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' dateofexp '\' analysisName];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' dateofexp]);
        mkdir(analysisName)
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(dateofexp,analysisName)
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

figNameBase = analysisName;

%%
frameratems = expt(iexp).frame_rate/1000;
cyctimems = cycTime/frameratems;
stimon_frames = input.nFramesOn;
stimoff_frames = input.nFramesOff;
stimontime = stimon_frames/frameratems;
stimofftime = stimoff_frames/frameratems;
trialresult = unique(trialOutcome);
responsedelay_frames = 3;

%% find sets of cells
DirFolder = expt(iexp).dirtuning;
run('cellSetsLG.m')
drivencells = find(any(mean(cycDataDFoverF_cmlvNoTarget{end},3) >0.05,1));% 

%% set cycle length
CYC = find(cycles == 5);
tempdata = cycDataDFoverF_cmlvNoTarget{CYC};
cells = [114 135 138 189];

%% plot figure for select cycle data, all cells, all trials, aud vs vis
figure;
for i = 1:length(cells)
    subplot(length(cells),1,i)
    icell = cells(i);
yminmax = [-0.02 0.1];
V_cycInd = cycV_ind{CYC};
AV_cycInd = cycAV_ind{CYC};
V_avg = mean(mean(tempdata(:,icell,V_cycInd),3),2);
AV_avg = mean(mean(tempdata(:,icell,AV_cycInd),3),2);
errbar_V = (std(mean(tempdata(:,icell,V_cycInd),2),[],3))/(sqrt(size(tempdata(:,icell,V_cycInd),3)));
errbar_AV = (std(mean(tempdata(:,icell,AV_cycInd),2),[],3))/(sqrt(size(tempdata(:,icell,AV_cycInd),3)));
xaxisMS = ([1:length(V_avg)]-prepress_frames)/frameratems;

% errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
p1 = shadedErrorBar(xaxisMS,AV_avg,errbar_AV,'k');
hold on

% errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
p2 = shadedErrorBar(xaxisMS,V_avg,errbar_V,'g');
hold on
% vline(0,'k')

hold on
startpatch = patch([0 stimontime stimontime 0],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
% set(startpatch,'FaceAlpha',0.15);
% set(startpatch,'EdgeColor','none');

% hold on
% offpatch = patch([stimontime stimofftime+stimontime stimofftime+stimontime stimontime],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
% set(offpatch,'FaceAlpha',0.5);
% set(offpatch,'EdgeColor','none');


hold on
for i = 1:cycles(CYC)-1
    L = (i*cyctimems);
    cyclepatch = patch([L L+stimontime L+stimontime L],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
    set(cyclepatch,'FaceAlpha',0.15);
    set(cyclepatch,'EdgeColor','none');
%     vline(L,'k:');
    hold on
end
vline((cycles(CYC)*cyctimems),'c');
hold on

ylim(yminmax)
text(-750, 0.008, {['cell# ' num2str(icell)]; ['vis trials = ' num2str(length(V_cycInd))]; ['aud trials = ' num2str(length(AV_cycInd))]})
xlabel('time (ms)')
ylabel('dF/F')
% axis square

end
suptitle({[mouse '-' dateofexp '-' runstr];[num2str(CYC) ' cycles']; ['mean driven cells, success trials']; 'gr = vis, bl = aud'})
% print([fnouttemp '\' mouse '-' dateofexp '_' num2str(CYC) ' cycles_selectcells'], '-dpdf');

%% plot tuned response to targets
run('FSAV_targetresp.m')
CYC = 5;
tempdata = cycDataDFoverF{CYC};

%find 90deg trials for CYC
Dir90 = find(DirectionDeg == Dirs(end));
Dir0 = find(DirectionDeg == Dirs(1));
CYCx = find(tCyclesOn == cycles(CYC));
CYCxDir90 = find(ismember(CYCx,Dir90));
CYCxDir0 = find(ismember(CYCx,Dir0)); 



%% mean response of cells to each target stim

targetMean = zeros(60,nCells,length(Dirs));
targetSTD = zeros(60,nCells,length(Dirs));
targetMeanResp = zeros(length(Dirs),nCells);
targetErrResp = zeros(length(Dirs),nCells);
for i = 1:length(Dirs)
    tempdatadir = dirTargetDataDFoverF{i};
    targetMean(:,:,i) = mean(tempdatadir,3);
    targetSTD(:,:,i) = std(tempdatadir,[],3);
    temptargetBL = squeeze(mean(targetMean(17:23,:,i),1));
    temptargetRSP = squeeze(mean(targetMean(24:27,:,i),1));
    targetMeanResp(i,:) = bsxfun(@minus,temptargetRSP,temptargetBL);
    targetErrResp(i,:) = bsxfun(@rdivide,squeeze(mean(targetSTD(24:27,:,i),1)),sqrt(size(tempdatadir,3)));
end

for i = 1:length(Dirs)
    ux = DirectionDeg==Dirs(i);
    dircounts(i) = sum(ux(strcmp(trialOutcome,'success')));
end

targetMean_xcells = zeros(length(Dirs),size(cellsSelect,2));
targetErr_xcells = zeros(length(Dirs),size(cellsSelect,2));
for i = 1:size(cellsSelect,2)
   cellspref = cellsSelect{i};
   targetMean_xcells(:,i) = mean(targetMeanResp(:,cellspref),2);
   targetErr_xcells(:,i) = squeeze(std(mean(targetSTD(24:27,cellspref,:),1),[],2)).\sqrt(dircounts)';
       
end

%% 
figure;
for i = 1:length(cells)
    subplot(length(cells)+1,1,i)
    icell = cells(i);
yminmax = [-0.03 0.35];
V_cycInd = CYCxDir90;
AV_cycInd = CYCxDir0;
V_avg = mean(mean(tempdata(:,icell,V_cycInd),3),2);
AV_avg = mean(mean(tempdata(:,icell,AV_cycInd),3),2);
errbar_V = (std(mean(tempdata(:,icell,V_cycInd),2),[],3))/(sqrt(size(tempdata(:,icell,V_cycInd),3)));
errbar_AV = (std(mean(tempdata(:,icell,AV_cycInd),2),[],3))/(sqrt(size(tempdata(:,icell,AV_cycInd),3)));
xaxisMS = ([1:length(V_avg)]-prepress_frames)/frameratems;

% errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
p1 = shadedErrorBar(xaxisMS,AV_avg,errbar_AV,'k');
hold on

% errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
p2 = shadedErrorBar(xaxisMS,V_avg,errbar_V,'g');
hold on
% vline(0,'k')

hold on
startpatch = patch([0 stimontime stimontime 0],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
% set(startpatch,'FaceAlpha',0.15);
% set(startpatch,'EdgeColor','none');

% hold on
% offpatch = patch([stimontime stimofftime+stimontime stimofftime+stimontime stimontime],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
% set(offpatch,'FaceAlpha',0.5);
% set(offpatch,'EdgeColor','none');


hold on
for i = 1:cycles(CYC)-1
    L = (i*cyctimems);
    cyclepatch = patch([L L+stimontime L+stimontime L],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
    set(cyclepatch,'FaceAlpha',0.15);
    set(cyclepatch,'EdgeColor','none');
%     vline(L,'k:');
    hold on
end
L = (cycles(CYC)*cyctimems);
targetpatch = patch([L L+stimontime L+stimontime L],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'c');
% vline((cycles(CYC)*cyctimems),'c');
hold on

ylim(yminmax)
text(-750, 0.008, {['cell# ' num2str(icell)]; ['vis trials = ' num2str(length(V_cycInd))]; ['aud trials = ' num2str(length(AV_cycInd))]})
xlabel('time (ms)')
ylabel('dF/F')
% axis square

end
suptitle({[mouse '-' dateofexp '-' runstr];[num2str(CYC) ' cycles']; ['mean driven cells, success trials']; 'gr = vis, bl = aud'})
print([fnouttemp '\' mouse '-' dateofexp '_' num2str(CYC) ' cycles_selectcells_withtarget'], '-dpdf');

%%
col = {'y','c'};
legendInfo = [];

figure;
for icell = 1:length(cells)
    errorbar(Dirs,targetMeanResp(:,cells(icell)),targetErrResp(:,cells(icell)),[col{icell} 'o'],'MarkerSize', 10)
    legendInfo{icell} = ['Cell ' num2str(cells(icell))];
    hold on
end
ylim([-0.1 0.4]) 
xlim([-10 Dirs(end)+5])
xlabel('Direction')
ylabel('dF/F')
set(gca,'XTick',Dirs)
legend(legendInfo)
title('Select Cell Tuning')
axis square
print([fnouttemp '\' mouse '-' dateofexp '_' num2str(CYC) ' cycles_selectcells_withtarget'], '-dpdf');

%% plot response to target averaged accross like-tuned cells.
col = {'k'; 'b'; 'r'; 'g'};
legendInfo = [];

figure;
for iori = 1:length(Oris)
%     errorbar(Dirs,targetMean_xcells(:,iori),targetErr_xcells(:,iori),[col{iori} 'o'],'MarkerSize', 10)
    plot(Dirs,targetMean_xcells(:,iori),[col{iori} 'o'],'MarkerSize', 10)
    legendInfo{iori} = [num2str(Oris(iori)) ' deg'];
    hold on
end
ylim([-0.01 0.05]) 
xlim([-10 Dirs(end)+5])
xlabel('Direction')
ylabel('dF/F')
set(gca,'XTick',Dirs)
legend(legendInfo)
title('Avg response to target accross like-tuned cells')
axis square
print([fnouttemp '\' mouse '-' dateofexp '_' num2str(CYC) ' cycles_avgresptotarget_sortedbycellpref'], '-dpdf');
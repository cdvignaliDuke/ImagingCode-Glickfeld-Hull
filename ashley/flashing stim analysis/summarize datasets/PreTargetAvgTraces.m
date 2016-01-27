iexp =1;
awFSAVdatasets_AL
cellgroupname = 'nCells';
figName = ['avgDFoverFpretarget_' cellgroupname '_success'];

ialign = 1;
run('divideupdatabyalignment.m')
%% find sets of cells
DirFolder = expt(iexp).dirtuning;
run('cellSets.m')
drivencells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetAvgTraces'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetAvgTraces')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetAvgTraces')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

%%
if length(eval(cellgroupname)) == 1
    cells = 1:eval(cellgroupname);
else
    cells = eval(cellgroupname);
end
%%
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
%     V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'failure')));
%     AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'failure')));
    V_avg = mean(mean(dataDFoverF(:,cells,V_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cells,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,cells,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,V_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,cells,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,AV_cycInd),3)));
    
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
    
    if icyc == 1
        title([num2str(length(cells)) ' cells'])
    else
    title({[num2str(length(V_cycInd)) ' visual trials; ']; [num2str(length(AV_cycInd)) ' vis+aud trials']})
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
%     ylim([-0.05 0.05])
end
suptitle([mouse '; ' date '; ' cellgroupname])

print([fnout ['\' figName '.pdf']], '-dpdf')
%%
if length(cells) > 12
    c = randperm(length(cells),12);
    cells = cells(c);
end

CYC = 6;

subpC = 4;
if length(cells) >= subpC;
    subpR = ceil(length(cells)/subpC);
else
    subpR = 1;
end

cellsavgFig = figure;
cellsalltrialsFig = figure;
for icell = 1:length(cells)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{CYC};
    V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
    AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
%     V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'ignore')));
%     AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'ignore')));
%     V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'failure')));
%     AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'failure')));
    V_avg = squeeze(mean(dataDFoverF(:,cells(icell),V_cycInd),3));
    AV_avg = squeeze(mean(dataDFoverF(:,cells(icell),AV_cycInd),3));
    errbar_V = (std(dataDFoverF(:,cells(icell),V_cycInd),[],3))/(sqrt(size(dataDFoverF(:,cells(icell),V_cycInd),3)));
    errbar_AV = (std(dataDFoverF(:,cells(icell),AV_cycInd),[],3))/(sqrt(size(dataDFoverF(:,cells(icell),AV_cycInd),3)));
    
    figure(cellsavgFig)
    subplot(subpR,subpC,icell);
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
    
    if icyc == 1
        title([num2str(length(cells)) ' cells'])
    else
    title({[num2str(length(V_cycInd)) ' visual trials; ']; [num2str(length(AV_cycInd)) ' vis+aud trials']; ['cell#' num2str(cells(icell))] })
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
%     ylim([-0.05 0.05])

    figure(cellsalltrialsFig)
    subplot(subpR,subpC,icell);
    plot(squeeze(dataDFoverF(:,cells(icell),V_cycInd)),'g')
    hold on
    plot(squeeze(dataDFoverF(:,cells(icell),AV_cycInd)),'k')
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
        title([num2str(length(cells)) ' cells'])
    else
    title({[num2str(length(V_cycInd)) ' visual trials; ']; [num2str(length(AV_cycInd)) ' vis+aud trials']; ['cell#' num2str(cells(icell))] })
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
end
suptitle([mouse '; ' date '; ' cellgroupname])

figure(cellsavgFig)
print([fnout ['\avgDFoverFpretarget_' cellgroupname '.pdf']], '-dpdf')

figure(cellsalltrialsFig)
print([fnout ['\allTrialspretarget_' cellgroupname '.pdf']], '-dpdf')

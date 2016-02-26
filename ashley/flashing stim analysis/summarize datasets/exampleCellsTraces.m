ms = 'AW14';
sn = '614';
t = '1239';
dt = '150623';
dirtuning = '006';
rettuning = '005';
imouse = 2;
iexp = 1;
cell_excAnt = 107;
cell_inhAnt = 161;
cell_excTar = 131;
cell_nr = 11;
cell_mat = [cell_excAnt cell_inhAnt cell_excTar cell_nr];
cellType_str = {'Exc - Ant';'Inh - Ant';'Exc - Tar';'NR'};
taskPart_str = {'Anticipation';'Target';'dir tuning, DG';'ret tuning, DG'};
pressAlign = 1;
catchAlign = 3;
targetAlign = 2;
cFA = 3;
cCR = 4;
hits = 1;
misses = 2;
n = 4;
n2 = 4;

close all
av = behavParamsAV;
dataGroup = ['awFSAVdatasets' datasetStr];
eval(dataGroup)
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
cycTime = mouse(1).expt(1).info(1).cyc_time;
exampleCellsFig = figure;

fnout = fullfile(rc.caOutputDir, dataGroup, [date '_' mouse_str]); %% maybe lose mouse_str


%% set params for figures
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
set(0,'DefaultaxesFontSize', 16)

tt = -pre_event_frames:post_event_frames-1;
baseStimFrames = 0:cycTime:post_event_frames-1;
baseStimFramesPreTar = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% anticipation period
figure(exampleCellsFig)
iCellSubplotXpos = [1 5 9 13];
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   plot(tt,squeeze(mean(mouse(imouse).expt(iexp).align(pressAlign).av(1).outcome(1).resp(:,cell_mat(icell),:),3)),'g');
   hold on
   plot(tt,squeeze(mean(mouse(imouse).expt(iexp).align(pressAlign).av(2).outcome(1).resp(:,cell_mat(icell),:),3)),'k');
   vline(baseStimFrames,':k')
   xlim([-10 post_event_frames])
   ylim([-0.15 0.15])
   ylabel({cellType_str{icell},'dF/F'})
   if icell == 1
       title(taskPart_str{1})
   end
end

%% target period
figure(exampleCellsFig)
iCellSubplotXpos = [2 6 10 14];
exptDirs = mouse(imouse).expt(iexp).visTargets;
colorsT = brewermap(length(exptDirs)+15,'YlGn');
colorindT = [3:2:length(exptDirs)+12];
colorsT = colorsT(colorindT(1:length(exptDirs)),:);
cellTarResp = [];
dirLegend = [];
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   for idir = 1:length(exptDirs)
       if idir == 1
           dirCol = [0 0 0];
       else
           dirCol = colorsT(idir,:);
       end
   cellTarResp(idir) = plot(tt,squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{idir}(:,cell_mat(icell),:),3)),'color',dirCol);
   dirLegend{idir} = num2str(exptDirs(idir));
   hold on
   end
%    plot(tt,squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(2).outcome(1).resp(:,cell_mat(icell),:),3)),'k');
   vline(baseStimFramesPreTar,':k')
   xlim([-10 20])
   ylim([-0.15 0.15])
%    ylabel(cellType_str{icell})
   if icell == 1
       title(taskPart_str{2})
       legend(cellTarResp,dirLegend,'Location','best')
   end
       vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
        hold on
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
        hold on
end

%% ori tuning curves
cellOriAndRetTuningCurves

figure(exampleCellsFig)
iCellSubplotXpos = [3 7 11 15];
for icell = 1:length(cell_mat)
    oriP = oriPref_ind(icell);
    dirP = dirPref_ind(icell);
    osi = chop(cellOSI(icell),2);
    dsi = chop(cellDSI(icell),2);
    subplot(n,n2,iCellSubplotXpos(icell))
    errorbar(directions,dFoverFDirResp(:,icell),errbarDirResp(:,icell),'ko-')
    hold on
    set(gca,'XTick',directions)
    xlabel('direction')
    xlim([-10 360])
    cellPref(1) = plot(directions(oriP),dFoverFDirResp(oriP,icell),'ro');
    hold on
    cellPref(2) = plot(directions(dirP),dFoverFDirResp(dirP,icell),'bo');
    legend(cellPref,{['OSI = ' num2str(osi)];['DSI = ' num2str(dsi)]},'location','northeast');
   if icell == 1
       title(taskPart_str{3})
   end
end

%% save
figure(exampleCellsFig)
print([fnout 'exampleCells' datasetStr '.pdf'], '-dpdf')
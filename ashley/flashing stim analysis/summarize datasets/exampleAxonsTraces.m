ms = 'AW16';
sn = '616';
t = '1056';
dt = '151130';
dirtuning = '007';
rettuning = '005';
imouse = 1;
iexp = 1;
% cell_excAnt = 24;
% cell_inhAnt = 14;
% cell_excTar = 332;
% cell_nr = 11;
cell_excAnt = mouse(imouse).expt(iexp).cells(1).ind(2);
cell_inhAnt = mouse(imouse).expt(iexp).cells(9).ind(1);
cell_excTar = mouse(imouse).expt(iexp).cells(13).ind(3);
cell_nr = 1;
cell_mat = [cell_excAnt cell_inhAnt cell_excTar cell_nr];
cellType_str = {'Exc - Ant';'Inh - Ant';'Exc - Tar';'NR'};
taskPart_str = {'Trial Start';'Anticipation';'Target';'dir tuning, DG';'ret tuning, DG'};
pressAlign = 1;
catchAlign = 3;
targetAlign = 2;
cFA = 3;
cCR = 4;
hits = 1;
misses = 2;
n = 4;
n2 = 4;

datasetStr = '_V1axonsPM';
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
cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;

fnout = fullfile(rc.caOutputDir, dataGroup, [date '_' mouse_str]); %% maybe lose mouse_str


%% set params for figures
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
set(0,'DefaultaxesFontSize', 8)

tt = -pre_event_frames:post_event_frames-1;
ttMs = tt/(cycTime/cycTimeMs);
baseStimFrames = 0:cycTime:post_event_frames-1;
baseStimFramesPreTar = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;
%%
exampleCellsFig = figure;
%% first response
figure(exampleCellsFig)
iCellSubplotXpos = [1 5 9 13];
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   tempV = mouse(imouse).expt(iexp).align(pressAlign).av(1).outcome(1).resp(:,cell_mat(icell),:);
   tempA = mouse(imouse).expt(iexp).align(pressAlign).av(2).outcome(1).resp(:,cell_mat(icell),:);
   mV = squeeze(mean(tempV,3));
   steV = std(tempV,[],3)/sqrt(size(tempV,3));
   mA = squeeze(mean(tempA,3));
   steA = std(tempA,[],3)/sqrt(size(tempA,3));
%    plot(ttMs,mV,'g');
   shadedErrorBar(ttMs,mV,steV,'g',0);
   hold on
%    plot(ttMs,mA,'k');
   shadedErrorBar(ttMs,mA,steA,'k',0);
   hold on
   xlim([-300 600])
   ylim([-0.05 0.16])
   vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
   hold on
   vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
   hold on
   ylabel({cellType_str{icell},'dF/F'})
   if icell == 1
       title(taskPart_str{1})
   end
   text(-250,0.1,['cell ' num2str(cell_mat(icell))]);
end

%% anticipation period
figure(exampleCellsFig)
iCellSubplotXpos = iCellSubplotXpos+1;
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   tempV = mouse(imouse).expt(iexp).align(pressAlign).av(1).outcome(1).resp(:,cell_mat(icell),:);
   tempA = mouse(imouse).expt(iexp).align(pressAlign).av(2).outcome(1).resp(:,cell_mat(icell),:);
   mV = squeeze(mean(tempV,3));
   steV = std(tempV,[],3)/sqrt(size(tempV,3));
   mA = squeeze(mean(tempA,3));
   steA = std(tempA,[],3)/sqrt(size(tempA,3));
   shadedErrorBar(ttMs,mV,steV,'g',0)
%    plot(ttMs,mV,'g');
   hold on
   shadedErrorBar(ttMs,mA,steA,'k',0)
%    plot(ttMs,mA,'k');
   xlim([-300 (post_event_frames/(cycTime/cycTimeMs))])
   ylim([-0.05 0.16])
   vline(baseStimFrames/(cycTime/cycTimeMs),':k')
   vline(0,'k')
   ylabel({cellType_str{icell},'dF/F'})
   if icell == 1
       title(taskPart_str{2})
   end
%    text(0,-0.1,['cell ' num2str(cell_mat(icell))]);
end

%% combine aud and vis trials
exCellsFig2 = figure;
% first response
figure(exCellsFig2)
iCellSubplotXpos = [1 5 9 13];
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   tempR = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(1).outcome(1).resp(:,cell_mat(icell),:),mouse(imouse).expt(iexp).align(pressAlign).av(2).outcome(1).resp(:,cell_mat(icell),:));
   mR = squeeze(mean(tempR,3));
   steR = std(tempR,[],3)/sqrt(size(tempR,3));
%    plot(ttMs,mA,'k');
   shadedErrorBar(ttMs,mR,steR,'k',0);
   hold on
   xlim([-300 600])
   ylim([-0.075 0.2])
   vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
   hold on
   vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
   hold on
   ylabel({cellType_str{icell},'dF/F'})
   if icell == 1
       title(taskPart_str{1})
   end
   text(-250,0.1,['cell ' num2str(cell_mat(icell))]);
end
% anticipation period
figure(exCellsFig2)
iCellSubplotXpos = iCellSubplotXpos+1;
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
   tempR = cat(3,mouse(imouse).expt(iexp).align(pressAlign).av(1).outcome(1).resp(:,cell_mat(icell),:),mouse(imouse).expt(iexp).align(pressAlign).av(2).outcome(1).resp(:,cell_mat(icell),:));
   mR = squeeze(mean(tempR,3));
   steR = std(tempR,[],3)/sqrt(size(tempR,3));
   shadedErrorBar(ttMs,mR,steR,'k',0);
   hold on
%    plot(ttMs,mA,'k');
   xlim([-300 (post_event_frames/(cycTime/cycTimeMs))])
   ylim([-0.075 0.2])
   vline(baseStimFrames/(cycTime/cycTimeMs),':k')
   vline(0,'k')
   ylabel({cellType_str{icell},'dF/F'})
   if icell == 1
       title(taskPart_str{2})
   end
%    text(0,-0.1,['cell ' num2str(cell_mat(icell))]);
end
%% target period
figure(exampleCellsFig)
iCellSubplotXpos = iCellSubplotXpos+1;
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
       tRespTC = squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{idir}(:,cell_mat(icell),:),3));
       preWinResp = mean(tRespTC(pre_win));
       tRespTC = tRespTC-preWinResp;
   cellTarResp(idir) = plot(ttMs,tRespTC,'color',dirCol);
   dirLegend{idir} = num2str(exptDirs(idir));
   hold on
   end

   dirLegend{idir} = num2str(exptDirs(idir));
   
   xlim([-300 600])
   ylim([-0.05 0.16])
%    plot(ttMs,squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(2).outcome(1).resp(:,cell_mat(icell),:),3)),'k');
   vline(baseStimFramesPreTar/(cycTime/cycTimeMs),':k')
%    ylabel(cellType_str{icell})
   if icell == 1
       title(taskPart_str{3})
       legend(cellTarResp,dirLegend,'Location','best')
   end
       vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
        hold on
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
        hold on
end

figure(exCellsFig2)
cellTarResp = [];
dirLegend = [];
for icell = 1:length(cell_mat)
   subplot(n,n2,iCellSubplotXpos(icell))
%    for idir = 1:length(exptDirs)
%        if idir == 1
%            dirCol = [0 0 0];
%        else
%            dirCol = colorsT(idir,:);
%        end
%        tRespTC = squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{idir}(:,cell_mat(icell),:),3));
%        preWinResp = mean(tRespTC(pre_win));
%        tRespTC = tRespTC-preWinResp;
%    cellTarResp(idir) = plot(ttMs,tRespTC,'color',dirCol);
%    dirLegend{idir} = num2str(exptDirs(idir));
%    hold on
%    end
   
   tRespTC = squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{length(exptDirs)}(:,cell_mat(icell),:),3));
   tRespTC_err = squeeze(std(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{length(exptDirs)}(:,cell_mat(icell),:),[],3))/sqrt(size(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{length(exptDirs)}(:,cell_mat(icell),:),3));
   preWinResp = mean(tRespTC(pre_win));
   tRespTC = tRespTC-preWinResp;
   shadedErrorBar(ttMs,tRespTC,tRespTC_err,'k',0);
   
   xlim([-300 600])
   ylim([-0.075 0.2])
%    plot(ttMs,squeeze(mean(mouse(imouse).expt(iexp).align(targetAlign).av(2).outcome(1).resp(:,cell_mat(icell),:),3)),'k');
   vline(baseStimFramesPreTar/(cycTime/cycTimeMs),':k')
%    ylabel(cellType_str{icell})
%    if icell == 1
%        title(taskPart_str{3})
%        legend(cellTarResp,dirLegend,'Location','best')
%    end
       vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
        hold on
        vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
        hold on
end


%% ori tuning curves
load(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'cellsSelect.mat'));
directions = [0 45 90 135 180 225 270 315];

figure(exampleCellsFig)
iCellSubplotXpos = iCellSubplotXpos+1

for icell = 1:length(cell_mat)
    cellPref = [];
    oriP = ori_ind_all(cell_mat(icell));
    dirP = max_dir_ind(cell_mat(icell));
    osi = chop(OSI(cell_mat(icell)),2);
    dsi = chop(DSI(cell_mat(icell)),2);
    subplot(n,n2,iCellSubplotXpos(icell))
    errorbar(directions,dFoverF_DirResp_avg_rect(:,cell_mat(icell)),dFoverF_DirResp_sem_rect(:,cell_mat(icell)),'ko-')
    hold on
    set(gca,'XTick',directions)
    xlabel('direction')
    xlim([-10 360])
    ylim([-0.05 0.5])
    if ~isnan(oriP) & ~isnan(osi)
    cellPref(1) = plot(directions(oriP),dFoverF_DirResp_avg_rect(oriP,cell_mat(icell)),'ro');
    hold on
    end
    if ~isnan(dirP) & ~isnan(dsi)
    cellPref(2) = plot(directions(dirP),dFoverF_DirResp_avg_rect(dirP,cell_mat(icell)),'bo');
    end
    legend(cellPref,{['OSI = ' num2str(osi)];['DSI = ' num2str(dsi)]},'location','northeast');
   if icell == 1
       title(taskPart_str{4})
   end
end

figure(exCellsFig2)
for icell = 1:length(cell_mat)
    cellPref = [];
    oriP = ori_ind_all(cell_mat(icell));
    dirP = max_dir_ind(cell_mat(icell));
    osi = chop(OSI(cell_mat(icell)),2);
    dsi = chop(DSI(cell_mat(icell)),2);
    subplot(n,n2,iCellSubplotXpos(icell))
    errorbar(directions,dFoverF_DirResp_avg_rect(:,cell_mat(icell)),dFoverF_DirResp_sem_rect(:,cell_mat(icell)),'ko-')
    hold on
    set(gca,'XTick',directions)
    xlabel('direction')
    xlim([-10 360])
    ylim([-0.05 0.5])
    if ~isnan(oriP) & ~isnan(osi)
    cellPref(1) = plot(directions(oriP),dFoverF_DirResp_avg_rect(oriP,cell_mat(icell)),'ro');
    hold on
    end
    if ~isnan(dirP) & ~isnan(dsi)
    cellPref(2) = plot(directions(dirP),dFoverF_DirResp_avg_rect(dirP,cell_mat(icell)),'bo');
    end
    legend(cellPref,{['OSI = ' num2str(osi)];['DSI = ' num2str(dsi)]},'location','northeast');
   if icell == 1
       title(taskPart_str{4})
   end
end

%% save
figure(exampleCellsFig)
print([fnout 'exampleCells' datasetStr '.pdf'], '-dpdf')
figure(exCellsFig2)
print([fnout 'exampleCells_acrossTrialTypes' datasetStr '.pdf'], '-dpdf')
% %% imaging FOV with cell masks
load(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'mask&TCDir.mat'));
maxDFoverF = readtiff(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'maxDFoverF.tif'));
% F = readtiff(fullfile(fnout,'Fimg.tif'));
sb_calib_x = 555.23/size(maxDFoverF,2); %um per pixel
sb_calib_y = 233.56/size(maxDFoverF,1);

umL50 = ceil(50/sb_calib_x);
umW5 = ceil(5/sb_calib_y);

crop_x_ind = 50:size(maxDFoverF,2)-50;
crop_y_ind = 50:size(maxDFoverF,1)-50;

maxDFoverF_crop = double(maxDFoverF(crop_y_ind,crop_x_ind));
% F_crop = double(F(crop_y_ind,crop_x_ind));

writetiff(maxDFoverF_crop,[fnout 'exampleCellsIMG' datasetStr]);
% writetiff(F_crop,[fnout 'exampleCellsIMG_F' datasetStr]);

sb_x_ind = 300:300+umL50;
sb_y_ind = 150:150+umW5;

img_sb = zeros(size(maxDFoverF_crop));
img_sb(sb_y_ind,sb_x_ind) = 1;
writetiff(img_sb,[fnout 'exampleCellsIMGsb'])
writetiff(maxDFoverF_crop,[fnout 'maxDFoverF_crop']);

img_sb = zeros(size(maxDFoverF));
img_sb(sb_y_ind,sb_x_ind) = 1;
writetiff(img_sb,[fnout 'notcropppedIMGsb'])


figure;
imagesc(img_sb)
colormap gray

cellMap = zeros(size(maxDFoverF,1),size(maxDFoverF,2),length(cell_mat));
for icell = 1:length(cell_mat)
    c = reshape(ismember(mask_boutons(:),cell_mat(icell)),size(mask_boutons));
    c(c == 1) = cell_mat(icell);
    cellMap(:,:,icell) = c;
end
cellMap_crop = cellMap(crop_y_ind,crop_x_ind,:);
writetiff(cellMap_crop,[fnout 'exampleCellsIMGcellmask' datasetStr]);


figure;
for icell = 1:length(cell_mat)
    figure;
    imagesc(cellMap{icell})
end












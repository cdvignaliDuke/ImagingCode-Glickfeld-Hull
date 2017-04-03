ms = 'AW72';
sn = '672';
dt = '170303';
dirtuning = '003';
imouse = 8;
iexp = 2;
cellType_str = {'Exc - Ant';'Exc - Tar'};
taskPart_str = {'Trial Start';'Anticipation';'Target';'dir tuning, DG';'ret tuning, DG'};
pressAlign = 1;
catchAlign = 3;
targetAlign = 2;
cFA = 3;
cCR = 4;
hits = 1;
misses = 2;
n = 3;
n2 = 4;
%%
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
cycTime = mouse(imouse).expt(iexp).info.cyc_time;
cycTimeMs = mouse(imouse).expt(iexp).info.cyc_time_ms;

fnout = fullfile(rc.caOutputDir, dataGroup, [date '_' ms '_' dt '_']); 
%%
d = mouse(imouse).expt(iexp);
%%
cell_excAnt = 51;%
cell_excTar = 7;%
cell_mat = [cell_excAnt cell_excTar];
%% set params for figures
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
set(0,'DefaultaxesFontSize', 8)

tt = -pre_event_frames:post_event_frames-1;
ttMs = tt/(cycTime/cycTimeMs);
baseStimFrames = 0:cycTime:post_event_frames-1;
baseStimFramesPreTar = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% combine aud and vis trials
exCellsFig2 = figure;
% first response
start_x_lim = [-250 700];
start_y_lim = [-0.05 0.1];
d_first = cat(3,d.align(pressAlign).av(1).outcome(1).cmlvCycResp{2},d.align(pressAlign).av(2).outcome(1).cmlvCycResp{2});
ttMs_first = ttMs(1:size(d_first,1));
figure(exCellsFig2)
iCellSubplotXpos = [1 5 9 13];
for icell = 1:length(cell_mat)
   sp1 = subplot(n,n2,iCellSubplotXpos(icell));
   mR = squeeze(mean(d_first(:,cell_mat(icell),:),3));
   steR = squeeze(ste(d_first(:,cell_mat(icell),:),3));

   shadedErrorBar(ttMs_first,mR,steR,'k',0);
   hold on
   figXAxis(sp1,'time(s)',start_x_lim)
   figYAxis(sp1,{cellType_str{icell},'dF/F'},start_y_lim)
   figAxForm(sp1)
   vline(0,'k--')
   if icell == 1
       title(taskPart_str{1})
   end
   text(-250,0.1,['cell ' num2str(cell_mat(icell))]);
end

% anticipation period
anti_x_lim = [-250 ((cycTime*8)/(cycTime/cycTimeMs))];
tr_y_lim = [-0.15 0.42];
d_anti = cat(3,d.align(pressAlign).av(1).outcome(1).cmlvCycResp{8},d.align(pressAlign).av(2).outcome(1).cmlvCycResp{8});
ttMs_anti = ttMs(1:size(d_anti,1));
figure(exCellsFig2)
iCellSubplotXpos = iCellSubplotXpos+1;
for icell = 1:length(cell_mat)
   sp2 = subplot(n,n2,iCellSubplotXpos(icell));

   mR = squeeze(mean(d_anti(:,cell_mat(icell),:),3));
   steR = squeeze(ste(d_anti(:,cell_mat(icell),:),3));
   shadedErrorBar(ttMs_anti,mR,steR,'k',0);
   hold on
   
   figXAxis(sp2,'time(s)',anti_x_lim);
   figYAxis(sp2,'dF/F',tr_y_lim);
   figAxForm(sp2);
   vline(0,'k--')

   if icell == 1
       title(taskPart_str{2})
   end
   
end

%target
tar_x_lim = [-250 500];
d_tar = mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{end};

iCellSubplotXpos = iCellSubplotXpos+1;

figure(exCellsFig2)
for icell = 1:length(cell_mat)
   sp3 = subplot(n,n2,iCellSubplotXpos(icell));
   tRespTC = squeeze(mean(d_tar(:,cell_mat(icell),:),3));
   tRespTC_err = squeeze(ste(d_tar(:,cell_mat(icell),:),3));

%    preWinResp = mean(tRespTC(pre_win));
%    tRespTC = tRespTC-preWinResp;
   shadedErrorBar(ttMs,tRespTC,tRespTC_err,'k',0);
   
   
   figXAxis(sp3,'time(s)',tar_x_lim);
   figYAxis(sp3,'dF/F',tr_y_lim);
   figAxForm(sp3);
   vline(0,'k--')

   if icell == 1
       title(taskPart_str{3})
   end
end


% ori tuning curves
load(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'cellsSelect.mat'));
directions = [0 45 90 135 180 225 270 315];

% figure(exampleCellsFig)
iCellSubplotXpos = iCellSubplotXpos+1;
ori_x_lim = [-10 360];
ori_y_lim = [-0.05 0.2];
figure(exCellsFig2)
for icell = 1:length(cell_mat)
    cellPref = [];
    oriP = ori_ind_all(cell_mat(icell));
    dirP = max_dir_ind(cell_mat(icell));
    osi = chop(OSI(cell_mat(icell)),2);
    dsi = chop(DSI(cell_mat(icell)),2);
    
    sp4 = subplot(n,n2,iCellSubplotXpos(icell));
    h = errorbar(directions,dFoverF_DirResp_avg_rect(:,cell_mat(icell)),dFoverF_DirResp_sem_rect(:,cell_mat(icell)),'ko-');
    h.MarkerFaceColor = [1 1 1];
    hold on
    
   figXAxis(sp4,'direction',ori_x_lim,directions);
   figYAxis(sp4,'dF/F',ori_y_lim);
   figAxForm(sp4);
    
%     if ~isnan(oriP) & ~isnan(osi)
%     cellPref(1) = plot(directions(oriP),dFoverF_DirResp_avg_rect(oriP,cell_mat(icell)),'ro');
%     hold on
%     end
%     if ~isnan(dirP) & ~isnan(dsi)
%     cellPref(2) = plot(directions(dirP),dFoverF_DirResp_avg_rect(dirP,cell_mat(icell)),'bo');
%     end
%     legend(cellPref,{['OSI = ' num2str(osi)];['DSI = ' num2str(dsi)]},'location','northeast');
   if icell == 1
       title(taskPart_str{4})
   end
end

%% response windows
oneS_fr = expt(ex).frame_rate;
oneStim_fr = oneS_fr/10;

start_stim = ones(size(ttMs_first));
start_stim(pre_event_frames+1:pre_event_frames+oneStim_fr) = 0;
start_img = repmat(start_stim,10,1);

start_lim_fr = zeros(1,2);
for i = 1:2
    [dum ind] = min(abs(ttMs-start_x_lim(i)));
    start_lim_fr(i) = ind;
end

anti_stim = ones(size(ttMs_anti));
anti_stim_ind = linspaceNDim(1:cycTime:cycTime*8,oneStim_fr:cycTime:cycTime*8,oneStim_fr);
anti_stim(anti_stim_ind(:)+pre_event_frames) = 0;
anti_img = repmat(anti_stim,10,1);

anti_lim_fr = zeros(1,2);
for i = 1:2
    [dum ind] = min(abs(ttMs-anti_x_lim(i)));
    anti_lim_fr(i) = ind;
end

tar_stim = ones(1,size(d_tar,1));
tar_stim(pre_event_frames+1:pre_event_frames+oneStim_fr) = 0;

tar_lim_fr = zeros(1,2);
for i = 1:2
    [dum ind] = min(abs(ttMs-tar_x_lim(i)));
    tar_lim_fr(i) = ind;
end

iCellSubplotXpos = [9 10 11];

figure(exCellsFig2)
colormap gray

sp = subplot(n,n2,9);
imagesc(start_img)
figXAxis(sp,'time (s)',start_lim_fr)
figAxForm(sp)

sp = subplot(n,n2,10);
imagesc(anti_stim)
figXAxis(sp,'time (s)',anti_lim_fr)
figAxForm(sp)

sp = subplot(n,n2,11);
imagesc(tar_stim)
figXAxis(sp,'time (s)',tar_lim_fr)
figAxForm(sp)

%% save
figure(exCellsFig2)
print([fnout 'exampleCells_acrossTrialTypes' datasetStr '.pdf'], '-dpdf', '-fillpage')
%% imaging FOV with cell masks
fnin = fullfile(rc.ashleyAnalysis,ms,expt(iexp).folder,dt);
load(fullfile(fnin,'max_images_crop.mat'))    
load(fullfile(fnin,'final_mask.mat')) 

maxDFoverF = max(cat(3,bx_crop,dir_crop),[],3);   

% % % % % 
% % % % % load(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'mask&TCDir.mat'));
% % % % % maxDFoverF = readtiff(fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,dirtuning,'maxDFoverF.tif'));
% % % % % % F = readtiff(fullfile(fn,'Fimg.tif'));
sb_calib_x = 555.23/size(maxDFoverF,2); %um per pixel
sb_calib_y = 233.56/size(maxDFoverF,1);

umL50 = ceil(50/sb_calib_x);
umW5 = ceil(5/sb_calib_y);

% % crop_x_ind = 50:size(maxDFoverF,2)-50;
% % crop_y_ind = 50:size(maxDFoverF,1)-50;
% % 
% % maxDFoverF_crop = double(maxDFoverF(crop_y_ind,crop_x_ind));
% % F_crop = double(F(crop_y_ind,crop_x_ind));

% % writetiff(maxDFoverF_crop,[fnout 'exampleCellsIMG' datasetStr]);
% % writetiff(F_crop,[fnout 'exampleCellsIMG_F' datasetStr]);

sb_x_ind = 500:500+umL50;
sb_y_ind = 220:220+umW5;

img_sb = zeros(size(maxDFoverF));
img_sb(sb_y_ind,sb_x_ind) = 1;
figure;colormap gray
imagesc(img_sb)

writetiff(img_sb,[fnout 'exampleCellsIMGsb'])
writetiff(maxDFoverF,[fnout,'exampleCellsIMG'])



cellMap = zeros(size(maxDFoverF,1),size(maxDFoverF,2),length(cell_mat));
for icell = 1:length(cell_mat)
    c = reshape(ismember(mask_cell(:),cell_mat(icell)),size(mask_cell));
    c(c == 1) = cell_mat(icell);
    cellMap(:,:,icell) = c;
end
% % cellMap_crop = cellMap(crop_y_ind,crop_x_ind,:);

writetiff(cellMap,[fnout 'exampleCellsIMGcellmask' datasetStr]);


figure;
for icell = 1:length(cell_mat)
    figure;
    imagesc(cellMap(:,:,icell))
end












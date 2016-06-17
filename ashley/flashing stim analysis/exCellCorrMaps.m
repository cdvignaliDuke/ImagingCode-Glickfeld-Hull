datasetStr = ['_V1'];
cellsInd = 14;

av = behavParamsAV;
eval(['awFSAVdatasets' datasetStr])
rc = behavConstsAV;
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end
if strcmp(rc.name,'ashley')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' mouse_str]); %% maybe lose mouse_str

pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
minTrialLengthFrames = mouse(1).expt(1).info.minTrialLengthFrames;
ialign = 1;
cycTime = mouse(1).expt(1).info(1).cyc_time;

ialign = 1;
% tt = -pre_event_frames:minTrialLengthFrames-1;
baseStimFrames = pre_event_frames:cycTime:minTrialLengthFrames+pre_event_frames-1;

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% example dataset 1
imouse = 1;
iexp = 2;
% for imouse = 1:3
cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;

% get responses and ns corr 
resp_int_V_b = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(1:pre_event_frames,cell_ind,:)))';
resp_int_V_1 = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(pre_event_frames:ceil(pre_event_frames+(minTrialLengthFrames/2)),cell_ind,:)))';
resp_int_V_2 = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(ceil(pre_event_frames+(minTrialLengthFrames/2)):pre_event_frames+minTrialLengthFrames,cell_ind,:)))';
nsc_v_b = corrcoef(resp_int_V_b);
nsc_v_1 = corrcoef(resp_int_V_1);
nsc_v_2 = corrcoef(resp_int_V_2);

resp_int_A_b = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(1:pre_event_frames,cell_ind,:)))';
resp_int_A_1 = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(pre_event_frames:ceil(pre_event_frames+(minTrialLengthFrames/2)),cell_ind,:)))';
resp_int_A_2 = squeeze(trapz(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(ceil(pre_event_frames+(minTrialLengthFrames/2)):pre_event_frames+minTrialLengthFrames,cell_ind,:)))';
nsc_a_b = corrcoef(resp_int_A_b);
nsc_a_1 = corrcoef(resp_int_A_1);
nsc_a_2 = corrcoef(resp_int_A_2);

%load cell mask
exptN = intersect(find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)),find(strcmp({expt.date},mouse(imouse).expt(iexp).date)));
load(fullfile(rc.ashleyAnalysis,expt(exptN).mouse,'two-photon imaging',expt(exptN).date,expt(exptN).dirtuning,'mask&TCDir.mat'))

nCells = mouse(imouse).expt(iexp).info.nCells;
rndCells = randperm(length(cell_ind),4);
[M,N] = size(mask_cell);
mask_vec = mask_cell(:);
mask_exclCells = ismember(mask_vec,setdiff(1:nCells,cell_ind));
mask_vec(mask_exclCells) = 0;

mask_vec_cellmat_v_b = repmat(mask_vec,1,4);
mask_vec_cellmat_v_1 = repmat(mask_vec,1,4);
mask_vec_cellmat_v_2 = repmat(mask_vec,1,4);
mask_vec_cellmat_a_b = repmat(mask_vec,1,4);
mask_vec_cellmat_a_1 = repmat(mask_vec,1,4);
mask_vec_cellmat_a_2 = repmat(mask_vec,1,4);
for i = 1:length(cell_ind)
    ind = find(mask_vec == cell_ind(i));
    for icell = 1:4
    val_v_b = nsc_v_b(i,rndCells(icell));
    val_v_1 = nsc_v_1(i,rndCells(icell));
    val_v_2 = nsc_v_2(i,rndCells(icell));
    val_a_b = nsc_a_b(i,rndCells(icell));
    val_a_1 = nsc_a_1(i,rndCells(icell));
    val_a_2 = nsc_a_2(i,rndCells(icell));
    mask_vec_cellmat_v_b(ind,icell) = val_v_b;
    mask_vec_cellmat_v_1(ind,icell) = val_v_1;
    mask_vec_cellmat_v_2(ind,icell) = val_v_2;
    mask_vec_cellmat_a_b(ind,icell) = val_a_b;
    mask_vec_cellmat_a_1(ind,icell) = val_a_1;
    mask_vec_cellmat_a_2(ind,icell) = val_a_2;
    end
end

mask_nscorr_4cells_v_b = reshape(mask_vec_cellmat_v_b,M,N,4);
mask_nscorr_4cells_v_1 = reshape(mask_vec_cellmat_v_1,M,N,4);
mask_nscorr_4cells_v_2 = reshape(mask_vec_cellmat_v_2,M,N,4);
mask_nscorr_4cells_a_b = reshape(mask_vec_cellmat_a_b,M,N,4);
mask_nscorr_4cells_a_1 = reshape(mask_vec_cellmat_a_1,M,N,4);
mask_nscorr_4cells_a_2 = reshape(mask_vec_cellmat_a_2,M,N,4);

mask_v_b = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date ' vis baseline']);
colormap(brewermap([],'*RdBu'))
mask_v_1 = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name ' vis 1']);
colormap(brewermap([],'*RdBu'))
mask_v_2 = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name ' vis 2']);
colormap(brewermap([],'*RdBu'))
mask_a_b = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date ' aud baseline']);
colormap(brewermap([],'*RdBu'))
mask_a_1 = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name ' aud 1']);
colormap(brewermap([],'*RdBu'))
mask_a_2 = figure;
suptitle([mouse(imouse).expt(iexp).mouse_name ' aud 2']);
colormap(brewermap([],'*RdBu'))
for iplot = 1:4
    figure(mask_v_b)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_v_b(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
    figure(mask_v_1)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_v_1(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
    figure(mask_v_2)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_v_2(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
    figure(mask_a_b)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_a_b(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
    figure(mask_a_1)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_a_1(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
    figure(mask_a_2)
    subplot(2,2,iplot)
    imagesc(mask_nscorr_4cells_a_2(:,:,iplot))
    title(['cell# ' num2cell(rndCells(iplot))])
    colorbar
    caxis([-1 1])
end
figure(mask_v_b)
print([fnout 'vis_nsCorrMap_exCell_baseline' datasetStr '.pdf'], '-dpdf')
figure(mask_a_b)
print([fnout 'aud_nsCorrMap_exCell_baseline' datasetStr '.pdf'], '-dpdf')

% scatter ns corr of cell pair on aud trials vs vis trials (2nd half only)
figure;
scatter(nsc_v_2(:),nsc_a_2(:),'k.')
xlim([-1 1]);
ylim([-1 1]);
hold on
plot([min(cat(1,nsc_v_2(:),nsc_a_2(:))):0.1:max(cat(1,nsc_v_2(:),nsc_a_2(:)))],[min(cat(1,nsc_v_2(:),nsc_a_2(:))):0.1:max(cat(1,nsc_v_2(:),nsc_a_2(:)))],'k--')
axis square


% view average timecourse of resp for aud and vis trials:
resp_tc_V = squeeze(mean(mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:),3))';
resp_tc_A = squeeze(mean(mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:),3))';

figure;
subplot(1,3,1)
colormap(brewermap([],'*RdBu'))
imagesc(resp_tc_V)
colorbar
axis square
caxis([-0.1 0.1])
vline(baseStimFrames,'k')
title({[mouse(imouse).expt(iexp).mouse_name '-' mouse(imouse).expt(iexp).date];'vis trial timecourse'})
subplot(1,3,2)
imagesc(resp_tc_A)
colorbar
axis square
caxis([-0.1 0.1])
vline(baseStimFrames,'k')
title('aud trial timecourse')
subplot(1,3,3)
imagesc(resp_tc_V-resp_tc_A)
colorbar
axis square
caxis([-0.1 0.1])
vline(baseStimFrames,'k')
title('vis-aud trial timecourse')


% overlay cell tcs that are highly correlated or anti-correlated
nplots = 16;
[hiCorrInd_r_a_2 hiCorrInd_c_a_2] = find((tril(nsc_a_2,-1) > 0.6));
[loCorrInd_r_a_2 loCorrInd_c_a_2] = find((tril(nsc_a_2,-1) < -0.6));
[hiCorrInd_r_a_b hiCorrInd_c_a_b] = find((tril(nsc_a_b,-1) > 0.6));
[loCorrInd_r_a_b loCorrInd_c_a_b] = find((tril(nsc_a_b,-1) < -0.6));

mostCorr_a_2 = mode(cat(1,hiCorrInd_r_a_2,hiCorrInd_c_a_2));
mostAntiCorr_a_2 = mode(cat(1,loCorrInd_r_a_2,loCorrInd_c_a_2));
mostCorr_a_b = mode(cat(1,hiCorrInd_r_a_b,hiCorrInd_c_a_b));
mostAntiCorr_a_b = mode(cat(1,loCorrInd_r_a_b,loCorrInd_c_a_b));

if length(hiCorrInd_r_a_2) > nplots
    hiCellsInd = randperm(length(hiCorrInd_r_a_2),nplots);
else
    hiCellsInd = 1:length(hiCorrInd_r_a_2);
end
if length(loCorrInd_r_a_2) > nplots
    loCellsInd = randperm(length(loCorrInd_r_a_2),nplots);
else
    loCellsInd = 1:length(loCorrInd_r_a_2);
end

if length(hiCorrInd_r_a_b) > nplots
    hiCellsInd_b = randperm(length(hiCorrInd_r_a_b),nplots);
else
    hiCellsInd_b = 1:length(hiCorrInd_r_a_b);
end
if length(loCorrInd_r_a_b) > nplots
    loCellsInd_b = randperm(length(loCorrInd_r_a_b),nplots);
else
    loCellsInd_b = 1:length(loCorrInd_r_a_b);
end

figure;
for iplot = 1:length(hiCellsInd)
    subplot(4,4,iplot)
    plot(resp_tc_A(hiCorrInd_r_a_2(hiCellsInd(iplot)),:),'r','linewidth',3)
    hold on
    plot(resp_tc_A(hiCorrInd_c_a_2(hiCellsInd(iplot)),:),'m','linewidth',3)
    title(['cells: ' num2str(hiCorrInd_r_a_2(hiCellsInd(iplot))) ',' num2str(hiCorrInd_c_a_2(hiCellsInd(iplot))) '; nsc_a_2: ' num2str(nsc_a_2(hiCorrInd_r_a_2(hiCellsInd(iplot)),hiCorrInd_c_a_2(hiCellsInd(iplot))))])
    if any(cat(1,hiCorrInd_r_a_2(hiCellsInd(iplot)),hiCorrInd_c_a_2(hiCellsInd(iplot))) == mostCorr_a_2)
        vline(100,'k')
    end
end
figure;
for iplot = 1:length(loCellsInd)
    subplot(4,4,iplot)
    plot(resp_tc_A(loCorrInd_r_a_2(loCellsInd(iplot)),:),'b','linewidth',3)
    hold on
    plot(resp_tc_A(loCorrInd_c_a_2(loCellsInd(iplot)),:),'c','linewidth',3)
    title(['cells: ' num2str(loCorrInd_r_a_2(loCellsInd(iplot))) ',' num2str(loCorrInd_c_a_2(loCellsInd(iplot))) '; nsc_a_2: ' num2str(nsc_a_2(loCorrInd_r_a_2(loCellsInd(iplot)),loCorrInd_c_a_2(loCellsInd(iplot))))])
    if any(cat(1,loCorrInd_r_a_2(loCellsInd(iplot)),loCorrInd_c_a_2(loCellsInd(iplot))) == mostAntiCorr_a_2)
        vline(100,'k')
    end
end

%scatter hi and anti correlated neurons responses (or residuals?)
resp_tc_Vall = mouse(imouse).expt(iexp).align(ialign).av(1).outcome(1).cmlvResp(:,cell_ind,:);
resp_tc_Aall = mouse(imouse).expt(iexp).align(ialign).av(2).outcome(1).cmlvResp(:,cell_ind,:);

resp_int_mean_A_2 = mean(resp_int_A_2,1);
resp_int_res_A_2 = bsxfun(@minus,resp_int_A_2,resp_int_mean_A_2);
resp_int_mean_A_b = mean(resp_int_A_b,1);
resp_int_res_A_b = bsxfun(@minus,resp_int_A_b,resp_int_mean_A_b);

figure;
for iplot = 1:length(hiCellsInd_b)
    subplot(4,4,iplot)
    scatter(resp_int_res_A_b(:,hiCorrInd_r_a_b(hiCellsInd_b(iplot))),resp_int_res_A_b(:,hiCorrInd_c_a_b(hiCellsInd_b(iplot))),'k.')
    xlim([-1 1]);
    ylim([-1 1]);
    hold on
    plot([min(min(resp_int_A_b)):0.1:max(max(resp_int_A_b))],[min(min(resp_int_A_b)):0.1:max(max(resp_int_A_b))],'k--')
    title(['cells: ' num2str(hiCorrInd_r_a_b(hiCellsInd_b(iplot))) ',' num2str(hiCorrInd_c_a_b(hiCellsInd_b(iplot))) '; nsc_a_b: ' num2str(nsc_a_b(hiCorrInd_r_a_b(hiCellsInd_b(iplot)),hiCorrInd_c_a_b(hiCellsInd_b(iplot))))])
    axis square
end
figure;
for iplot = 1:length(loCellsInd_b)
    subplot(4,4,iplot)
    scatter(resp_int_res_A_b(:,loCorrInd_r_a_b(loCellsInd_b(iplot))),resp_int_res_A_b(:,loCorrInd_c_a_b(loCellsInd_b(iplot))),'k.')
    xlim([-1 1]);
    ylim([-1 1]);
    hold on
    plot([min(min(resp_int_A_b)):0.1:max(max(resp_int_A_b))],[min(min(resp_int_A_b)):0.1:max(max(resp_int_A_b))],'k--')
    title(['cells: ' num2str(loCorrInd_r_a_b(loCellsInd_b(iplot))) ',' num2str(loCorrInd_c_a_b(loCellsInd_b(iplot))) '; nsc_a_b: ' num2str(nsc_a_b(loCorrInd_r_a_b(loCellsInd_b(iplot)),loCorrInd_c_a_b(loCellsInd_b(iplot))))])
    axis square
end
%% example dataset 2
imouse = 2;
iexp = 1;

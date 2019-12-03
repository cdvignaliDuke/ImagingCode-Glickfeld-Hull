clear all
close all
ds = 'szTuning_PM';
eval(ds)
%%
calib = 1/26.6; %mm per pixel

rc = behavConstsAV;
% imgParams_FSAV
fnout = 'Z:\home\ashley\Manuscripts\Size Tuning\Matlab Analysis';

%% analysis params

nbaselinefr = round((params.nBaselineMs/1000)*params.frameRate);
tt_stimOnS_label = -0.75:0.25:0.75;
tcStartTimeS = -1;

nexp = size(expt,2);
%% load and align eye tc for each experiment

eyeExpt = struct;
eyeExpt_ret = struct;
for iexp = 1:nexp
    
%     load eye data
    ms = expt(iexp).mouse;
    dt = expt(iexp).date;
    szFolder = expt(iexp).sizeTuningFolder{1};
    szTime = expt(iexp).sizeTuningTime{1};
    fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,szFolder);
    
    load(fullfile(fn,'eyeTC.mat'))
    
    eyeExpt(iexp).mouse = ms;
    eyeExpt(iexp).date = dt;
    
%     load and extract stim info
    mw = loadMworksFile(ms,dt,szTime);
    
    on = mw.nScansOn;
    off = mw.nScansOff;
    sz = celleqel2mat_padded(mw.tGratingDiameterDeg);
    if iexp == 1
        tt_stimOnFr = (1:(off+on))-(off+params.nFramesVisDelay_VSR);
        tt_stimOnS = double(tt_stimOnFr)./params.frameRate;
        respwin = (off+1):(off+on);
    end
    
    basewin = (off-nbaselinefr+1):off;
    
    nfr = length(Area)-1;
    ntrials = floor(nfr./(on+off));
    if ntrials > length(sz)
        ntrials = length(sz);
    end
    nfr_tr = (on+off)*ntrials;
    Area = Area(1:nfr_tr);
    Centroid = Centroid(1:nfr_tr,:);
    
    ind = ~isnan(Area);
    Area_interp = interp1(find(ind),Area(ind),1:length(Area));
    Centroid_interp = nan(size(Centroid));
    Centroid_interp(:,1) = interp1(find(ind),Centroid(ind,1),1:length(Area));
    Centroid_interp(:,2) = interp1(find(ind),Centroid(ind,2),1:length(Area));
    
    area_tr = reshape(Area_interp,[on+off,ntrials]).*...
        params.eyeCalibMmPerPix.*1000./params.eyeCalibUmPerDeg;
    centr_tr = reshape(Centroid_interp,[on+off,ntrials,2]).*...
        params.eyeCalibMmPerPix.*1000./params.eyeCalibUmPerDeg;
    
    eyeExpt(iexp).area.allTrialsTC = area_tr./mean(area_tr(basewin,:),1);
        
    eyeExpt(iexp).pos(1).name = 'horizontal';
    eyeExpt(iexp).pos(2).name = 'vertical';
    eyeExpt(iexp).pos(1).allTrialsTC_noBL = centr_tr(:,:,1);
    eyeExpt(iexp).pos(2).allTrialsTC_noBL = centr_tr(:,:,2);
    eyeExpt(iexp).pos(1).allTrialsTC = centr_tr(:,:,1) - ...
        mean(squeeze(centr_tr(basewin,:,1)),1);
    eyeExpt(iexp).pos(2).allTrialsTC = centr_tr(:,:,2) - ...
        mean(squeeze(centr_tr(basewin,:,2)),1);
    
    d = squeeze(mean(centr_tr(respwin,:,:)));
    nt = size(d,1);
    distFromMeanPupil = nan(nt,1);
    for itrial = 1:nt
        distFromMeanPupil(itrial) = pdist([d(itrial,:);mean(d,1)]);
    end
    
    eyeExpt(iexp).pos(3).name = 'relative distance';
    eyeExpt(iexp).pos(3).position = distFromMeanPupil;
        
    eyeExpt(iexp).tSz = sz;
    
    
%     load eye data from retinotopy
    eyeExpt_ret(iexp).mouse = ms;
    eyeExpt_ret(iexp).date = dt;
    retFolder = expt(iexp).retinotopyFolder{1};
    retTime = expt(iexp).retinotopyTime{1};
    fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,retFolder);
    mw = loadMworksFile(ms,dt,retTime);
    az = celleqel2mat_padded(mw.tGratingDiameterDeg);
    el = celleqel2mat_padded(mw.tGratingDiameterDeg);
    
    if exist(fullfile(fn,'eyeTC.mat'),'file') == 2
        load(fullfile(fn,'eyeTC.mat'))  
    
    %     load and extract stim info

        on = mw.nScansOn;
        off = mw.nScansOff;

        basewin = (off-nbaselinefr+1):off;

        nfr = length(Area)-1;
        ntrials = floor(nfr./(on+off));
        if ntrials > length(az)
            ntrials = length(az);
        end
        nfr_tr = (on+off)*ntrials;
        Area = Area(1:nfr_tr);
        Centroid = Centroid(1:nfr_tr,:);

        ind = ~isnan(Area);
        Area_interp = interp1(find(ind),Area(ind),1:length(Area));
        Centroid_interp = nan(size(Centroid));
        Centroid_interp(:,1) = interp1(find(ind),Centroid(ind,1),1:length(Area));
        Centroid_interp(:,2) = interp1(find(ind),Centroid(ind,2),1:length(Area));

        area_tr = reshape(Area_interp,[on+off,ntrials]).*...
            params.eyeCalibMmPerPix.*1000./params.eyeCalibUmPerDeg;
        centr_tr = reshape(Centroid_interp,[on+off,ntrials,2]).*...
            params.eyeCalibMmPerPix.*1000./params.eyeCalibUmPerDeg;
        eyeExpt_ret(iexp).area.allTrialsTC = area_tr./mean(area_tr(basewin,:),1);        
        eyeExpt_ret(iexp).pos(1).name = 'horizontal';
        eyeExpt_ret(iexp).pos(2).name = 'vertical';
        eyeExpt_ret(iexp).pos(1).allTrialsTC_noBL = centr_tr(:,:,1);
        eyeExpt_ret(iexp).pos(2).allTrialsTC_noBL = centr_tr(:,:,2);
        eyeExpt_ret(iexp).pos(1).allTrialsTC = centr_tr(:,:,1) - ...
            mean(squeeze(centr_tr(basewin,:,1)),1);
        eyeExpt_ret(iexp).pos(2).allTrialsTC = centr_tr(:,:,2) - ...
            mean(squeeze(centr_tr(basewin,:,2)),1);

        respwin_ret = (off+1):(off+on);
        d = squeeze(mean(centr_tr(respwin_ret,:,:)));
        nt = size(d,1);
        distFromMeanPupil = nan(nt,1);
        for itrial = 1:nt
            distFromMeanPupil(itrial) = pdist([d(itrial,:);mean(d,1)]);
        end

        eyeExpt_ret(iexp).pos(3).position = distFromMeanPupil;

    else
        ntrials = length(az);
        eyeExpt_ret(iexp).pos(3).name = 'relative distance';
        eyeExpt_ret(iexp).pos(3).position = nan(ntrials,1);
        
    end
    eyeExpt_ret(iexp).tAz = az;
    eyeExpt_ret(iexp).tEl = el;
    
end

relDist_stdDev = nan(1,nexp);
for iexp = 1:nexp
    relDist_stdDev(iexp) = std(eyeExpt(iexp).pos(3).position);
end
fprintf('Standard Dev. of Relative Distance Moved: %s+/-%s\n',...
    num2str(round(mean(relDist_stdDev),2,'significant')),...
    num2str(round(std(relDist_stdDev),2,'significant')))

save(fullfile(fnout,'eyeSummary'),'eyeExpt','eyeExpt_ret')

%%

%% save example eye image and scalbar
exExpt = 2;

ms = expt(exExpt).mouse;
dt = expt(exExpt).date;
szFolder = expt(exExpt).sizeTuningFolder{1};
fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,szFolder);
d = readtiff(fullfile(fn,'eye_movie.tif'));
exFrame = randsample(1:size(d,3),1);
[c,r,metric] = imfindcircles(d(:,:,exFrame),expt(exExpt).eyeradrange);
[~,ind] = max(metric);
                

figure; colormap gray
subplot 121
imagesc(d(:,:,exFrame))
hold on
viscircles(c(ind,:),r(ind),'EdgeColor','w');
title([ms '-' dt '; frame ' num2str(exFrame)])
figAxForm

subplot 122
sb_img = zeros(size(d(:,:,exFrame)));
[ypix,xpix] = size(sb_img);
xpix_1mm = round(1/params.eyeCalibMmPerPix);
ypix_100um = round(0.1/params.eyeCalibMmPerPix);
sb_xInd = (1:xpix_1mm)+10;
sb_yInd = (1:ypix_100um)+10;
sb_img(sb_yInd,sb_xInd) = 1;
imagesc(sb_img)
title('Scalebar 1.0x0.1mm')
figAxForm

print(fullfile(fnout,'exEyeImage'),'-dpdf','-fillpage')
%% plot each mouse
areaLim = [0.6 1.4];
posLim = [-3 1];
distLim = [0 6];
distLimHist = [0 20];
distBinEdges = 0:15;
tcStartFr = find(tt_stimOnS > tcStartTimeS,1);

eyeSum = struct;
eyeSum.area.tc = cell(1,nexp);
eyeSum.area.szTC = cell(1,nexp);
eyeSum.pos(1).szTC = cell(1,nexp);
eyeSum.pos(2).szTC = cell(1,nexp);
eyeSum.pos(3).relPosAllTrials = cell(1,nexp);
eyeSum.pos(3).szRelPos = cell(1,nexp);
for iexp = 1:nexp
    tSz = round(eyeExpt(iexp).tSz,2,'significant');
    sizes = unique(tSz);
    nsz = length(sizes);
    
    areaSzTC = cell(1,nsz);
    hrzSzTC = cell(1,nsz);
    vrtSzTC = cell(1,nsz);
    hrzSzTC_noBL = cell(1,nsz);
    vrtSzTC_noBL = cell(1,nsz);
    posSz = cell(1,nsz);
    for isz = 1:nsz
        ind = tSz == sizes(isz);
        areaSzTC{isz} = eyeExpt(iexp).area.allTrialsTC(:,ind);
        hrzSzTC{isz} = eyeExpt(iexp).pos(1).allTrialsTC(:,ind);
        vrtSzTC{isz} = eyeExpt(iexp).pos(2).allTrialsTC(:,ind);
        hrzSzTC_noBL{isz} = eyeExpt(iexp).pos(1).allTrialsTC_noBL(:,ind);
        vrtSzTC_noBL{isz} = eyeExpt(iexp).pos(2).allTrialsTC_noBL(:,ind);
        posSz{isz} = eyeExpt(iexp).pos(3).position(ind);
    end
    
    eyeSum.pos(3).relPosAllTrials{iexp} = eyeExpt(iexp).pos(3).position;
    eyeSum.pos(3).szRelPos{iexp} = cellfun(@mean,posSz);
    
    figure
    suptitle([eyeExpt(iexp).mouse '-' eyeExpt(iexp).date])
    subplot 221
    x = tt_stimOnS(tcStartFr:end);
    y = mean(eyeExpt(iexp).area.allTrialsTC(tcStartFr:end,:),2);
    yerr = ste(eyeExpt(iexp).area.allTrialsTC(tcStartFr:end,:),2);
    h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'% Change', areaLim)
    figAxForm
    hold on
    hline(1,'k--')
    vline([respwin(1) respwin(end)]-tcStartFr,'r:')
    title('Pupil Area')
    eyeSum.area.tc{iexp} = y;
    
    subplot 223
    x = tt_stimOnS(tcStartFr:end);
    colors = parula(nsz+1);
    for isz = 1:nsz
        hold on
        y = mean(areaSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'% Change', areaLim)
    figAxForm
    hline(1,'k--')
    eyeSum.area.szTC{iexp} = cell2mat(cellfun(@(x) mean(x(tcStartFr:end,:),2),...
        areaSzTC,'unif',0));
    
    subplot 222
    x = sizes;
    for isz = 1:nsz
        hold on
        y = mean(mean(areaSzTC{isz}(respwin,:),2),1);
        yerr = ste(mean(areaSzTC{isz}(respwin,:),1),2);
        h = errorbar(x(isz),y,yerr,'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'% Change', areaLim)
    figAxForm
    hline(1,'k--')
    
    fn = fullfile(rc.ashleyAnalysis,eyeExpt(iexp).mouse,'two-photon imaging',...
        eyeExpt(iexp).date);
    print(fullfile(fn,'eyeSummary_area'),'-dpdf','-fillpage')
    
    figure 
    subplot 331
    hrzMean = mean(mean(eyeExpt(iexp).pos(1).allTrialsTC_noBL(respwin,:)));
    vrtMean = mean(mean(eyeExpt(iexp).pos(2).allTrialsTC_noBL(respwin,:)));
    for isz = 1:nsz
        hold on
        x = mean(hrzSzTC_noBL{isz}(respwin,:),1) - hrzMean;
        y = mean(vrtSzTC_noBL{isz}(respwin,:),1) - vrtMean;
        h = plot(x,y,'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    posLim_noBL = rangeMinMax(...
        cat(2,mean(eyeExpt(iexp).pos(1).allTrialsTC_noBL(respwin,:),1) - hrzMean,...
        mean(eyeExpt(iexp).pos(2).allTrialsTC_noBL(respwin,:),1) - vrtMean));
    figXAxis([],'Rel. Horiz. Pos. (deg)',posLim_noBL)
    figYAxis([],'Rel. Vert. Pos. (deg)',posLim_noBL)
    figAxForm
    vline(0,'k--')
    hline(0,'k--')
    
    subplot 334
    y = mean(eyeExpt(iexp).pos(1).allTrialsTC_noBL(respwin,:),1) - hrzMean;
    histogram(y,posLim_noBL(1):posLim_noBL(end))
    hold on
    figXAxis([],'Rel. Horiz. Pos. (deg)',posLim_noBL)
    figYAxis([],'N Trials',[])
    figAxForm
    vline(0,'k--')
    
    subplot 337
    y = mean(eyeExpt(iexp).pos(2).allTrialsTC_noBL(respwin,:),1) - vrtMean;
    histogram(y,posLim_noBL(1):posLim_noBL(end))
    hold on
    figXAxis([],'Rel. Vert. Pos. (deg)',posLim_noBL)
    figYAxis([],'N Trials',[])
    figAxForm
    vline(0,'k--')
    
    subplot 332 
    y = eyeExpt(iexp).pos(3).position;
    h = histogram(y,distBinEdges);
    figXAxis([],'Dist from mean pos (deg)',distLimHist)
    figYAxis([],'N Trials', [])
    figAxForm
    hold on
    vline(params.eyePosThreshold_deg,'r:')
    vline(mean(y),'k:')
    title(sprintf('Mean dist: %s',num2str(round(mean(y),2,'significant'))))
    
    subplot 333
    y = cellfun(@mean,posSz);
    yerr = cellfun(@(x) ste(x,1),posSz);
    for isz = 1:nsz
        hold on
        h = errorbar(sizes(isz),y(isz),yerr(isz),'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'Dist from mean pos (deg)', distLim)
    figAxForm
    
    subplot 335
    x = tt_stimOnS(tcStartFr:end);
    for isz = 1:nsz
        hold on
        y = mean(hrzSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Diff from Baseline(deg)', posLim)
    figAxForm
    hline(0,'k--')
    title('Horiz. Pupil Pos.')
    eyeSum.pos(1).szTC{iexp} = cell2mat(cellfun(@(x) mean(x(tcStartFr:end,:),2),...
        hrzSzTC,'unif',0));
    
    subplot 338
    x = sizes;
    for isz = 1:nsz
        hold on
        y = mean(mean(hrzSzTC{isz}(respwin,:),2),1);
        yerr = ste(mean(hrzSzTC{isz}(respwin,:),1),2);
        h = errorbar(x(isz),y,yerr,'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'Diff from Baseline(deg)', posLim)
    figAxForm
    hline(0,'k--')
    
    subplot 336
    x = tt_stimOnS(tcStartFr:end);
    for isz = 1:nsz
        hold on
        y = mean(vrtSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Diff from Baseline(deg)', posLim)
    figAxForm
    hline(0,'k--')
    title('Vert. Pupil Pos.')
    eyeSum.pos(2).szTC{iexp} = cell2mat(cellfun(@(x) mean(x(tcStartFr:end,:),2),...
        vrtSzTC,'unif',0));
    
    subplot 339
    x = sizes;
    for isz = 1:nsz
        hold on
        y = mean(mean(vrtSzTC{isz}(respwin,:),2),1);
        yerr = ste(mean(vrtSzTC{isz}(tcStartFr:end,:),1),2);
        h = errorbar(x(isz),y,yerr,'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'Diff from Baseline(deg)', posLim)
    figAxForm
    hline(0,'k--')
    
    fn = fullfile(rc.ashleyAnalysis,eyeExpt(iexp).mouse,'two-photon imaging',...
        eyeExpt(iexp).date);
    print(fullfile(fn,'eyeSummary_position'),'-dpdf','-fillpage')
end
%% plot across mice
figure
suptitle('All Expt')
subplot 331
x = tt_stimOnS(tcStartFr:end);
for iexp = 1:nexp
    hold on
    y = eyeSum.area.tc{iexp};
    plot(x,y,'-');
end
y = mean(cell2mat(eyeSum.area.tc),2);
yerr = ste(cell2mat(eyeSum.area.tc),2);
h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'% Change', areaLim)
figAxForm
hold on
hline(1,'k--')
vline([respwin(1) respwin(end)]-tcStartFr,'r:')
title('Pupil Area')
subplot 334
x = tt_stimOnS(tcStartFr:end);
colors = parula(nsz+1);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),eyeSum.area.szTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'% Change', areaLim)
figAxForm
hline(1,'k--')
subplot 337
x = sizes;
for isz = 1:nsz
    hold on
    y_all = mean(cell2mat(cellfun(@(x) x(respwin-tcStartFr,isz),eyeSum.area.szTC,'unif',0)),1);
    yerr = ste(y_all,2);
    h = plot(ones(1,nexp).*x(isz),y_all,'.','MarkerSize',10);
    h.Color = colors(isz,:);
    h = errorbar(x(isz),mean(y_all),yerr,'.','MarkerSize',10);
    h.Color = colors(isz,:);
end
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'% Change', areaLim)
figAxForm
hline(1,'k--')

subplot 332
for iexp = 1:nexp
    hold on
    y = eyeSum.pos(3).relPosAllTrials{iexp};
    h = histogram(y,distBinEdges);
end
figXAxis([],'Dist from mean pos (deg)',distLimHist)
figYAxis([],'N Trials', [])
figAxForm
hold on
vline(params.eyePosThreshold_deg,'r:')
y = cellfun(@mean,eyeSum.pos(3).relPosAllTrials);
vline(mean(y),'k:')
title(sprintf('Mean dist: %s+/-%s',num2str(round(mean(y),2,'significant')),...
    num2str(round(ste(y,2),2,'significant'))))

subplot 333
for iexp = 1:nexp
    y = eyeSum.pos(3).szRelPos{iexp};
%     for isz = 1:nsz
        hold on
        h = plot(sizes,y,'-');
        h.Color = [.5 .5 .5];
%     end
end
y = cell2mat(eyeSum.pos(3).szRelPos');
yerr = ste(y,1);
p=anova1(y,[],'off');
errorbar(sizes,mean(y),yerr,'k.-','MarkerSize',10)
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Dist from mean pos (deg)', distLim)
figAxForm
title(sprintf('%s+/-%s p=%s',num2str(round(mean(y(:)),2,'significant')),...
    num2str(round(ste(y(:),1),2,'significant')),...
    num2str(round(p,2,'significant'))))

subplot 335
x = tt_stimOnS(tcStartFr:end);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),eyeSum.pos(1).szTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Diff from Baseline(deg)', posLim)
figAxForm
hline(0,'k--')
title('Horiz. Pupil Pos.')
subplot 338
x = sizes;
for isz = 1:nsz
    hold on
    y_all = mean(cell2mat(cellfun(@(x) x(respwin-tcStartFr,isz),eyeSum.pos(1).szTC,'unif',0)),1);
    yerr = ste(y_all,2);
    h = plot(ones(1,nexp).*x(isz),y_all,'.','MarkerSize',10);
    h.Color = colors(isz,:);
    h = errorbar(x(isz),mean(y_all),yerr,'.','MarkerSize',10);
    h.Color = colors(isz,:);
end
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Diff from Baseline(deg)', posLim)
figAxForm
hline(0,'k--')

subplot 336
x = tt_stimOnS(tcStartFr:end);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),eyeSum.pos(2).szTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Diff from Baseline(deg)', posLim)
figAxForm
hline(0,'k--')
title('Vert. Pupil Pos.')
subplot 339
x = sizes;
for isz = 1:nsz
    hold on
    y_all = mean(cell2mat(cellfun(@(x) x(respwin-tcStartFr,isz),eyeSum.pos(2).szTC,'unif',0)),1);
    yerr = ste(y_all,2);
    h = plot(ones(1,nexp).*x(isz),y_all,'.','MarkerSize',10);
    h.Color = colors(isz,:);
    h = errorbar(x(isz),mean(y_all),yerr,'.','MarkerSize',10);
    h.Color = colors(isz,:);
end
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Diff from Baseline(deg)', posLim)
figAxForm
hline(0,'k--')
print(fullfile(fnout,'eyeSummary'),'-dpdf','-fillpage')
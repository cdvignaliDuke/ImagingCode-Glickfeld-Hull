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
    
    area_tr = reshape(Area_interp,[on+off,ntrials]).*params.eyeCalibMmPerPix;
    centr_tr = reshape(Centroid_interp,[on+off,ntrials,2]).*params.eyeCalibMmPerPix;
    
    eyeExpt(iexp).area.allTrialsTC = area_tr./mean(area_tr(basewin,:),1);
        
    eyeExpt(iexp).pos(1).name = 'horizontal';
    eyeExpt(iexp).pos(2).name = 'vertical';
    eyeExpt(iexp).pos(1).allTrialsTC = centr_tr(:,:,1) - ...
        mean(squeeze(centr_tr(basewin,:,1)),1);
    eyeExpt(iexp).pos(2).allTrialsTC = centr_tr(:,:,2) - ...
        mean(squeeze(centr_tr(basewin,:,2)),1);
    
    eyeExpt(iexp).tSz = sz;
    if iexp == 1
        tt_stimOnFr = (1:(off+on))-(off+params.nFramesVisDelay_VSR);
        tt_stimOnS = double(tt_stimOnFr)./params.frameRate;
        respwin = (off+4):(off+12);
    end
end

%% plot each mouse
areaLim = [0.6 1.4];
posLim = [-0.05 0.05];
tcStartFr = find(tt_stimOnS > tcStartTimeS,1);

eyeSum = struct;
eyeSum.area.tc = cell(1,nexp);
eyeSum.area.szTC = cell(1,nexp);
eyeSum.pos(1).tc = cell(1,nexp);
eyeSum.pos(1).szTC = cell(1,nexp);
eyeSum.pos(2).tc = cell(1,nexp);
eyeSum.pos(2).szTC = cell(1,nexp);
for iexp = 1:nexp
    tSz = round(eyeExpt(iexp).tSz,2,'significant');
    sizes = unique(tSz);
    nsz = length(sizes);
    
    areaSzTC = cell(1,nsz);
    hrzSzTC = cell(1,nsz);
    vrtSzTC = cell(1,nsz);
    for isz = 1:nsz
        ind = tSz == sizes(isz);
        areaSzTC{isz} = eyeExpt(iexp).area.allTrialsTC(:,ind);
        hrzSzTC{isz} = eyeExpt(iexp).pos(1).allTrialsTC(:,ind);
        vrtSzTC{isz} = eyeExpt(iexp).pos(2).allTrialsTC(:,ind);
    end
    
    
    msFig = figure;
    suptitle([eyeExpt(iexp).mouse '-' eyeExpt(iexp).date])
    subplot 331
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
    
    subplot 334
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
    
    subplot 337
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
    
    figure(msFig)
    subplot 332
    x = tt_stimOnS(tcStartFr:end);
    y = mean(eyeExpt(iexp).pos(1).allTrialsTC(tcStartFr:end,:),2);
    yerr = ste(eyeExpt(iexp).pos(1).allTrialsTC(tcStartFr:end,:),2);
    h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Abs. Change (mm)', posLim)
    figAxForm
    hold on
    hline(0,'k--')
    vline([respwin(1) respwin(end)]-tcStartFr,'r:')
    title('Horiz. Pupil Pos.')
    eyeSum.pos(1).tc{iexp} = y;
    
    subplot 335
    x = tt_stimOnS(tcStartFr:end);
    for isz = 1:nsz
        hold on
        y = mean(hrzSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Abs. Change (mm)', posLim)
    figAxForm
    hline(0,'k--')
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
    figYAxis([],'Abs. Change (mm)', posLim)
    figAxForm
    hline(0,'k--')
    
    
    figure(msFig)
    subplot 333
    x = tt_stimOnS(tcStartFr:end);
    y = mean(eyeExpt(iexp).pos(2).allTrialsTC(tcStartFr:end,:),2);
    yerr = ste(eyeExpt(iexp).pos(2).allTrialsTC(tcStartFr:end,:),2);
    h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Abs. Change', posLim)
    figAxForm
    hold on
    hline(0,'k--')
    vline([respwin(1) respwin(end)]-tcStartFr,'r:')
    title('Vert. Pupil Pos.')
    eyeSum.pos(2).tc{iexp} = y;
    
    subplot 336
    x = tt_stimOnS(tcStartFr:end);
    for isz = 1:nsz
        hold on
        y = mean(vrtSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Abs. Change (mm)', posLim)
    figAxForm
    hline(0,'k--')
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
    figYAxis([],'Abs. Change (mm)', posLim)
    figAxForm
    hline(0,'k--')
    
    fn = fullfile(rc.ashleyAnalysis,eyeExpt(iexp).mouse,'two-photon imaging',...
        eyeExpt(iexp).date);
    print(fullfile(fn,'eyeSummary'),'-dpdf','-fillpage')
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
x = tt_stimOnS(tcStartFr:end);
for iexp = 1:nexp
    hold on
    y = eyeSum.pos(1).tc{iexp};
    plot(x,y,'-');
end
y = mean(cell2mat(eyeSum.pos(1).tc),2);
yerr = ste(cell2mat(eyeSum.pos(1).tc),2);
h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hold on
hline(0,'k--')
vline([respwin(1) respwin(end)]-tcStartFr,'r:')
title('Horiz. Pupil Pos.')
subplot 335
x = tt_stimOnS(tcStartFr:end);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),eyeSum.pos(1).szTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hline(0,'k--')
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
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hline(0,'k--')

subplot 333
x = tt_stimOnS(tcStartFr:end);
for iexp = 1:nexp
    hold on
    y = eyeSum.pos(2).tc{iexp};
    plot(x,y,'-');
end
y = mean(cell2mat(eyeSum.pos(2).tc),2);
yerr = ste(cell2mat(eyeSum.pos(2).tc),2);
h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hold on
hline(0,'k--')
vline([respwin(1) respwin(end)]-tcStartFr,'r:')
title('Vert. Pupil Pos.')
subplot 336
x = tt_stimOnS(tcStartFr:end);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),eyeSum.pos(2).szTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hline(0,'k--')
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
figYAxis([],'Abs. Change (mm)', posLim)
figAxForm
hline(0,'k--')
print(fullfile(fnout,'eyeSummary'),'-dpdf','-fillpage')
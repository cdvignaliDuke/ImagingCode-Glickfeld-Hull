clear all
close all
ds = 'szTuning_PM';
eval(ds)
%%
rc = behavConstsAV;
% imgParams_FSAV
fnout = 'Z:\home\ashley\Manuscripts\Size Tuning\Matlab Analysis';

%% analysis params
runningThreshold_cps = 0;


nbaselinefr = round((params.nBaselineMs/1000)*params.frameRate);
tt_stimOnS_label = -0.75:0.25:0.75;
tcStartTimeS = -1;

nexp = size(expt,2);
%% load and align run tc for each experiment, new PM expts and all axon expts

runExpt = struct;
for iexp = 1:nexp
    
%     load mworks data
    ms = expt(iexp).mouse;
    dt = expt(iexp).date;
    szFolder = expt(iexp).sizeTuningFolder{1};
    szTime = expt(iexp).sizeTuningTime{1};
    fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,szFolder);
    
    mw = loadMworksFile(ms,dt,szTime);
    
    runExpt(iexp).mouse = ms;
    runExpt(iexp).date = dt;
    
%     calculate running speed on each frame
    ws_cps = wheelSpeedCalc(mw,32,'red');
    
    on = mw.nScansOn;
    off = mw.nScansOff;
    sz = celleqel2mat_padded(mw.tGratingDiameterDeg);
    
    basewin = (off-nbaselinefr+1):off;
    
    nfr = length(ws_cps)-1;
    ntrials = floor(nfr./(on+off));
    if ntrials > length(sz)
        ntrials = length(sz);
    end
    nfr_tr = (on+off)*ntrials;
    
    running = ws_cps(1:nfr_tr);
    
    running_tr = reshape(running,[on+off,ntrials]);
    
    runExpt(iexp).allTrialsTC = running_tr;
     
    runExpt(iexp).tSz = sz;
    if iexp == 1
        tt_stimOnFr = (1:(off+on))-(off+params.nFramesVisDelay_VSR);
        tt_stimOnS = double(tt_stimOnFr)./params.frameRate;
        respwin = (off+4):(off+12);
    end
end

area_list = {'LM', 'AL', 'PM'};
exptInd = [];
runExpt_axons = struct;
for iarea = 1:size(area_list,2)
    ds = ['szTuning_axons_' area_list{iarea}];
    eval(ds)
    nexp = size(expt,2);
    wheelOn = nan(1,nexp);
    for iexp = 1:nexp
        mouse = expt(iexp).mouse;
        expDate = expt(iexp).date;
        retRun = expt(iexp).retinotopyFolder;
        retTime = expt(iexp).retinotopyTime;
        szRun = expt(iexp).sizeTuningFolder;
        szTime = expt(iexp).sizeTuningTime{1};
        
        mwSize = loadMworksFile(mouse,expDate,szTime);
        if sum(cell2mat(mwSize.wheelSpeedValues))>0
            wheelOn(1,iexp) = true;
        else
            wheelOn(1,iexp) = false;
        end
        
        if wheelOn(iexp)
            if sum(wheelOn(~isnan(wheelOn))) == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
            ws_cps = wheelSpeedCalc(mwSize,32,'red');
            
            on = mwSize.nScansOn;
            off = mwSize.nScansOff;
            sz = celleqel2mat_padded(mwSize.tGratingDiameterDeg);

            basewin = (off-nbaselinefr+1):off;

            nfr = length(ws_cps)-1;
            ntrials = floor(nfr./(on+off));
            if ntrials > length(sz)
                ntrials = length(sz);
            end
            nfr_tr = (on+off)*ntrials;

            running = ws_cps(1:nfr_tr);
            running_tr = reshape(running,[on+off,ntrials]);
            
            runExpt_axons(exptN).allTrialsTC = running_tr;
            runExpt_axons(exptN).tSz = sz;
        else
        end        
    end
    exptInd = cat(2,exptInd,wheelOn);
end
  

save(fullfile(fnout,'runInfoEaExpt'),'runExpt')
%% plot each mouse
wsLim = [-2 2];
runProbLim = [0 0.15];
tcStartFr = find(tt_stimOnS > tcStartTimeS,1);

runSum = struct;
runSum.wheelSpeedTC = cell(1,nexp);
runSum.wheelSpeedSzTC = cell(1,nexp);
runSum.runProbSz = cell(1,nexp);
for iexp = 1:nexp
    tSz = round(runExpt(iexp).tSz,2,'significant');
    sizes = unique(tSz);
    nsz = length(sizes);
    
    runSzTC = cell(1,nsz);
    for isz = 1:nsz
        ind = tSz == sizes(isz);
        runSzTC{isz} = runExpt(iexp).allTrialsTC(:,ind);
    end
    runProbSz = cellfun(@(x) ...
        sum(sum(x(respwin,:)>runningThreshold_cps)>0)./length(x),...
        runSzTC);
    runSum.runProbSz{iexp} = runProbSz;
        
    msFig = figure;
    suptitle([runExpt(iexp).mouse '-' runExpt(iexp).date])
    subplot 221
    x = tt_stimOnS(tcStartFr:end);
    y = mean(runExpt(iexp).allTrialsTC(tcStartFr:end,:),2);
    yerr = ste(runExpt(iexp).allTrialsTC(tcStartFr:end,:),2);
    h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Wheel Speed (cm/s)', wsLim)
    figAxForm
    hold on
    hline(0,'k--')
    vline(tt_stimOnS([respwin(1) respwin(end)]),'r:')
    title('All Trials')
    runSum.wheelSpeedTC{iexp} = y;
    
    subplot 222
    x = tt_stimOnS(tcStartFr:end);
    colors = parula(nsz+1);
    for isz = 1:nsz
        hold on
        y = mean(runSzTC{isz}(tcStartFr:end,:),2);
        h = plot(x,y,'-');
        h.Color = colors(isz,:);
    end
    figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
    figYAxis([],'Wheel Speed (cm/s)', wsLim)
    figAxForm
    hline(0,'k--')
    runSum.wheelSpeedSzTC{iexp} = cell2mat(cellfun(@(x) mean(x(tcStartFr:end,:),2),...
        runSzTC,'unif',0));
    
    subplot 223
    x = sizes;
    for isz = 1:nsz
        hold on
        y = mean(mean(runSzTC{isz}(respwin,:),2),1);
        yerr = ste(mean(runSzTC{isz}(respwin,:),1),2);
        h = errorbar(x(isz),y,yerr,'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'Wheel Speed (cm/s)', wsLim)
    figAxForm
    hline(0,'k--')
    
    subplot 224
    x = sizes;
    for isz = 1:nsz
        hold on
        h = plot(sizes(isz),runProbSz(isz),'.','MarkerSize',10);
        h.Color = colors(isz,:);
    end
    figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
    figYAxis([],'Prob Run During Stim', runProbLim)
    figAxForm
    
    print(fullfile(fn,'runSummary'),'-dpdf','-fillpage')
end
%% plot across mice
figure
suptitle('All Expt')
subplot 221
x = tt_stimOnS(tcStartFr:end);
for iexp = 1:nexp
    hold on
    y = runSum.wheelSpeedTC{iexp};
    plot(x,y,'-');
end
y = mean(cell2mat(runSum.wheelSpeedTC),2);
yerr = ste(cell2mat(runSum.wheelSpeedTC),2);
h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Wheel Speed (cm/s)', wsLim)
figAxForm
hold on
hline(0,'k--')
vline(tt_stimOnS([respwin(1) respwin(end)]),'r:')
title('All Trials')

subplot 222
x = tt_stimOnS(tcStartFr:end);
colors = parula(nsz+1);
for isz = 1:nsz
    hold on
    y = mean(cell2mat(cellfun(@(x) x(:,isz),runSum.wheelSpeedSzTC,'unif',0)),2);
    h = plot(x,y,'-');
    h.Color = colors(isz,:);
end
figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
figYAxis([],'Wheel Speed (cm/s)', wsLim)
figAxForm
hold on
hline(0,'k--')
vline(tt_stimOnS([respwin(1) respwin(end)]),'r:')
title('All Trials')

subplot 223
x = sizes;
for isz = 1:nsz
    hold on
    y_all = mean(cell2mat(cellfun(@(x) x(respwin-tcStartFr,isz),runSum.wheelSpeedSzTC,'unif',0)),1);
    yerr = ste(y_all,2);
    h = plot(ones(1,nexp).*x(isz),y_all,'.','MarkerSize',10);
    h.Color = colors(isz,:);
    h = errorbar(x(isz),mean(y_all),yerr,'.','MarkerSize',10);
    h.Color = colors(isz,:);
end
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Wheel Speed (cm/s)', wsLim)
figAxForm
hline(0,'k--')

subplot 224
x = sizes;
for iexp = 1:nexp
    hold on
    h = plot(sizes,runSum.runProbSz{iexp},'-');
    h.Color = [0.5 0.5 0.5];
end
for isz = 1:nsz
    hold on
    y = cellfun(@(x) x(isz), runSum.runProbSz);
    h = plot(sizes(isz).*ones(1,nexp), ...
        y,'.','MarkerSize',10);
    h.Color = colors(isz,:);
    h = errorbar(sizes(isz),mean(y),ste(y,2),'.','MarkerSize',10);
    h.Color = colors(isz,:);
end
errorbar(sizes,mean(cell2mat(runSum.runProbSz')),...
    ste(cell2mat(runSum.runProbSz'),1),'k-');
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Prob Run During Stim', runProbLim)
figAxForm

print(fullfile(fnout,'runSummary'),'-dpdf','-fillpage')
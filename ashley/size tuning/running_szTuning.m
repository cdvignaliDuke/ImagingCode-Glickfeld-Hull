clear all
close all
ds = 'szTuning_PM';
eval(ds)
szParams = params;
%%
rc = behavConstsAV;
% imgParams_FSAV
fnout = 'Z:\home\ashley\Manuscripts\Size Tuning\Matlab Analysis';

%% analysis params
runningThreshold_cps = szParams.runSpeedThreshold_cps;


nbaselinefr = round((szParams.nBaselineMs/1000)*szParams.frameRate);
tt_stimOnS_label = -0.75:0.25:0.75;
tcStartTimeS = -1;

nexp = size(expt,2);
%% load and align run tc for each experiment

runExpt_cells = struct;
for iexp = 1:nexp
    
%     load mworks data
    ms = expt(iexp).mouse;
    dt = expt(iexp).date;
%     szFolder = expt(iexp).sizeTuningFolder{1};
    szTime = expt(iexp).sizeTuningTime{1};
%     fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,szFolder);
    
    mw = loadMworksFile(ms,dt,szTime);
    
    runExpt_cells(iexp).mouse = ms;
    runExpt_cells(iexp).date = dt;
    
%     calculate running speed on each frame
    ws_cps = wheelSpeedCalc(mw,32,'red');
    
    on = mw.nScansOn;
    off = mw.nScansOff;
    sz = celleqel2mat_padded(mw.tGratingDiameterDeg);
    if iexp == 1
        tt_stimOnFr = (1:(off+on))-(off+szParams.nFramesVisDelay_VSR);
        tt_stimOnS = double(tt_stimOnFr)./szParams.frameRate;
        respwin_sz = (off+1):(off+on);
    end
    
    basewin = (off-nbaselinefr+1):off;
    
    nfr = length(ws_cps)-1;
    ntrials = floor(nfr./(on+off));
    if ntrials > length(sz)
        ntrials = length(sz);
    end
    nfr_tr = (on+off)*ntrials;
    
    running = ws_cps(1:nfr_tr);
    
    running_tr = reshape(running,[on+off,ntrials]);
    
    runExpt_cells(iexp).sz.allTrialsTC = running_tr;
    runExpt_cells(iexp).sz.allTrialsSpeedCPS = mean(running_tr(respwin_sz,:),1);     
    runExpt_cells(iexp).sz.tSz = sz;
    
%     runFolder = expt(iexp).sizeTuningFolder{1};
    retTime = expt(iexp).retinotopyTime{1};
%     fn = fullfile(rc.ashleyAnalysis,ms,'two-photon imaging',dt,runFolder);    
    mw = loadMworksFile(ms,dt,retTime);
    
    ws_cps = wheelSpeedCalc(mw,32,'red');
    
    on = mw.nScansOn;
    off = mw.nScansOff;
    az = celleqel2mat_padded(mw.tGratingAzimuthDeg);
    el = celleqel2mat_padded(mw.tGratingAzimuthDeg);
    
    nfr = length(ws_cps)-1;
    ntrials = floor(nfr./(on+off));
    if ntrials > length(az)
        ntrials = length(az);
    end
    nfr_tr = (on+off)*ntrials;
    
    running = ws_cps(1:nfr_tr);
    
    running_tr = reshape(running,[on+off,ntrials]);
    respwin_ret = (off+1):(off+on);
    runExpt_cells(iexp).ret.allTrialsTC = running_tr;
    runExpt_cells(iexp).ret.allTrialsSpeedCPS = mean(running_tr(respwin_ret,:),1);     
    runExpt_cells(iexp).ret.tAz = az;
    runExpt_cells(iexp).ret.tEl = el;
end

area_list = {'LM', 'AL', 'PM'};
exptInd = [];
runExpt_temp_axons = struct;
for iarea = 1:size(area_list,2)
    clear expt
    ds = ['szTuning_axons_' area_list{iarea}];
    eval(ds)
    nexp = size(expt,2);
    wheelOn = nan(1,nexp);
    mouseOn = nan(1,nexp);
    for iexp = 1:nexp
        mouse = expt(iexp).mouse;
        if any(strcmp(mice2use,mouse))
            mouseOn(iexp) = true;
        else
            mouseOn(iexp) = false;
        end
        expDate = expt(iexp).date;
        retRun = expt(iexp).retinotopyFolder;
        retTime = expt(iexp).retinotopyTime{1};
        szRun = expt(iexp).sizeTuningFolder;
        szTime = expt(iexp).sizeTuningTime{1};
        
        mwSize = loadMworksFile(mouse,expDate,szTime);
        mwRet = loadMworksFile(mouse,expDate,retTime);
        if sum(abs(cell2mat(mwSize.wheelSpeedValues)))>0
            wheelOn(1,iexp) = true;
        else
            wheelOn(1,iexp) = false;
        end
        
        if iarea == 1 & iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        if wheelOn(iexp)
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
            
            runExpt_temp_axons(exptN).mouse = mouse;
            runExpt_temp_axons(exptN).date = expDate;
            runExpt_temp_axons(exptN).sz.allTrialsTC = running_tr;
            runExpt_temp_axons(exptN).sz.allTrialsSpeedCPS = mean(running_tr(respwin_sz,:),1);
            runExpt_temp_axons(exptN).sz.tSz = sz;
            
            ws_cps = wheelSpeedCalc(mwRet,32,'red');            
            on = mwRet.nScansOn;
            off = mwRet.nScansOff;
            az = celleqel2mat_padded(mwRet.tGratingAzimuthDeg);
            el = celleqel2mat_padded(mwRet.tGratingAzimuthDeg);
            nfr = length(ws_cps)-1;
            ntrials = floor(nfr./(on+off));
            if ntrials > length(az)
                ntrials = length(az);
            end
            nfr_tr = (on+off)*ntrials;

            running = ws_cps(1:nfr_tr);
            running_tr = reshape(running,[on+off,ntrials]);
            
            runExpt_temp_axons(exptN).ret.allTrialsTC = running_tr;
            runExpt_temp_axons(exptN).ret.allTrialsSpeedCPS = mean(running_tr(respwin_ret,:),1);
            runExpt_temp_axons(exptN).ret.tAz = az;
            runExpt_temp_axons(exptN).ret.tEl = el;
        else
            runExpt_temp_axons(exptN).mouse = mouse;
            runExpt_temp_axons(exptN).date = expDate;
        end        
    end
    exptInd = cat(2,exptInd,wheelOn & mouseOn);
end
runExpt_axons = runExpt_temp_axons(exptInd==1);
save(fullfile(fnout,'runSummary'),'runExpt_cells','runExpt_axons','szParams')
%% plot each mouse
doMsPlot = false;
wsLim = [-2 12];
runProbLim = [0 0.5];
tcStartFr = find(tt_stimOnS > tcStartTimeS,1);

runExpt = cat(2,runExpt_cells,runExpt_axons);
nexp = size(runExpt,2);
runSum = struct;
runSum.wheelSpeedTC = cell(1,nexp);
runSum.wheelSpeedSzTC = cell(1,nexp);
runSum.runProbSz = cell(1,nexp);
runSum.runProb = nan(1,nexp);
for iexp = 1:nexp
    fn = fullfile(rc.ashleyAnalysis,runExpt(iexp).mouse,...
        'two-photon imaging',runExpt(iexp).date);
    tSz = round(runExpt(iexp).sz.tSz,2,'significant');
    sizes = unique(tSz);
    nsz = length(sizes);
    
    runSzTC = cell(1,nsz);
    for isz = 1:nsz
        ind = tSz == sizes(isz);
        runSzTC{isz} = runExpt(iexp).sz.allTrialsTC(:,ind);
    end
    runProbSz = cellfun(@(x) ...
        sum(mean(x(respwin_sz,:),1)>runningThreshold_cps)./size(x,2),runSzTC);
    runSum.runProbSz{iexp} = runProbSz;
    runSum.runProb(iexp) = ...
        sum(mean(runExpt(iexp).sz.allTrialsTC(respwin_sz,:),1)>runningThreshold_cps)./...
        size(runExpt(iexp).sz.allTrialsTC,2);
        
%     if doMsPlot
        msFig = figure;
        suptitle([runExpt(iexp).mouse '-' runExpt(iexp).date])
        subplot 221
        x = tt_stimOnS(tcStartFr:end);
        y = mean(runExpt(iexp).sz.allTrialsTC(tcStartFr:end,:),2);
        yerr = ste(runExpt(iexp).sz.allTrialsTC(tcStartFr:end,:),2);
        h = shadedErrorBar_chooseColor(x,y,yerr,[0 0 0]);
        figXAxis([],'Time from Stim On (s)',[x(1) x(end)],tt_stimOnS_label,tt_stimOnS_label)
        figYAxis([],'Wheel Speed (cm/s)', wsLim)
        figAxForm
        hold on
        hline(0,'k--')
        vline(tt_stimOnS([respwin_sz(1) respwin_sz(end)]),'r:')
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
            y = mean(mean(runSzTC{isz}(respwin_sz,:),2),1);
            yerr = ste(mean(runSzTC{isz}(respwin_sz,:),1),2);
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
%     end
end

fprintf('Run Prob. Range: %s-%s\n',num2str(min(runSum.runProb)),...
    num2str(max(runSum.runProb)))
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
vline(tt_stimOnS([respwin_sz(1) respwin_sz(end)]),'r:')
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
vline(tt_stimOnS([respwin_sz(1) respwin_sz(end)]),'r:')
title('All Trials')

subplot 223
x = sizes;
for isz = 1:nsz
    hold on
    y_all = mean(cell2mat(cellfun(@(x) x(respwin_sz-tcStartFr,isz),runSum.wheelSpeedSzTC,'unif',0)),1);
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
% for isz = 1:nsz
%     hold on
%     y = cellfun(@(x) x(isz), runSum.runProbSz);
%     h = plot(sizes(isz).*ones(1,nexp), ...
%         y,'.','MarkerSize',10);
% %     h.Color = colors(isz,:);
% %     h = errorbar(sizes(isz),mean(y),ste(y,2),'.','MarkerSize',10);
% %     h.Color = colors(isz,:);
% end
errorbar(sizes,mean(cell2mat(runSum.runProbSz')),...
    ste(cell2mat(runSum.runProbSz'),1),'k-');
[p,~,stats] = anova1(cell2mat(runSum.runProbSz'),[],'off');
figXAxis([],'Stim. Size (deg)',[0 100],sizes,sizes)
figYAxis([],'Prob Run During Stim', runProbLim)
figAxForm
title(sprintf('Run Prob all trials: %s+/-%s, p=%s',...
    num2str(round(mean(runSum.runProb),2,'significant')),...
    num2str(round(ste(runSum.runProb,2),2,'significant')),...
    num2str(round(p,2,'significant'))))

print(fullfile(fnout,'runSummary'),'-dpdf','-fillpage')
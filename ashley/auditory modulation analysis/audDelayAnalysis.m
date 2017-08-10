clear all
close all
ds = '_V13trialtypes';
%%
rc = behavConstsAV;
eval(['awData_audMod' ds])

%% load data
fn = fullfile(rc.ashleyAnalysis,'Expt Summaries',['awData_audMod' ds]);
load(fullfile(fn,['awData_audMod' ds '_CaSummary']),'ms');

%% get ori tuning and fits
doFits = 1;
if doFits
    neuronTuning = getOriTuningAllCells(ms,fn);
    neuronTimeCourse = getTrialTimeCourseAllCells(ms,fn);
else
    load(fullfile(fn,'oriTuningAndFits.mat'))
    load(fullfile(fn,'neuronTrialTimeCourses.mat'))
end
%% some params
theta_cutoff = 30;
respWinS = .75;
preFrames = ms(1).pre_event_frames;
ntrtype = length(ms(1).delayTrials);
trType_str = ms(1).audDelayType;
frRateHz = ms(1).frRateHz;
% wheelRateHz = double(ms(2).wheelRateHz);
% preWheel = double(ms(2).pre_event_wheel);
oris = [];
expt_str = cell(1,size(ms,2));
for iexp = 1:size(ms,2)
    msOris = ms(iexp).orientations;
    oris = cat(2,oris,msOris);
    expt_str{iexp} = [ms(iexp).subnum '-' ms(iexp).date];
end   
oris = unique(oris);
nori = length(oris);
visOnly_ind = find(strcmp(trType_str,'vis only'));

nboot = 10;
  
pre_win = preFrames-round(respWinS*frRateHz)+1:preFrames;
resp_win = preFrames + round(.1*frRateHz):preFrames + round((respWinS+.1)*frRateHz);

c = brewermap(ntrtype+1,'*Blues');
delay_colors(1:ntrtype-1,:) = c(2:end-1,:);
delay_colors(ntrtype,:) = [0 0 0];
expt_colors = brewermap(size(ms,2),'*Dark2');

%% pupil
nexp = size(ms,2);
pexp = false(1,nexp);
pexpt_str = cell(1,nexp);
psz = cell(1,nexp);
horiz = cell(1,nexp);
vert = cell(1,nexp);
pszLong = cell(1,nexp);
horizLong = cell(1,nexp);
vertLong = cell(1,nexp);
for iexp = 1:nexp
    if ~isempty(ms(iexp).pupil{1}{1})
        p = ms(iexp).pupil{1};
        h = ms(iexp).pupil{2};
        v = ms(iexp).pupil{3};
        
        psz{iexp} = cellfun(@(x) nanmean(x,2),p,'unif',0);
        horiz{iexp} = cellfun(@(x) nanmean(x,2),h,'unif',0);
        vert{iexp} = cellfun(@(x) nanmean(x,2),v,'unif',0);
        
        p = ms(iexp).pupilLong{1};
        h = ms(iexp).pupilLong{2};
        v = ms(iexp).pupilLong{3};
        
        pszLong{iexp} = cellfun(@(x) nanmean(x,2),p,'unif',0);
        horizLong{iexp} = cellfun(@(x) nanmean(x,2),h,'unif',0);
        vertLong{iexp} = cellfun(@(x) nanmean(x,2),v,'unif',0);
        
        pexp(iexp) = true;
        expt_str{iexp} = [ms(iexp).subnum '-' ms(iexp).date];
    end    
end

p = psz(pexp);
hz = horiz(pexp);
vt = vert(pexp);

pupilLong = pszLong(pexp);
horizLong = horizLong(pexp);
vertLong = vertLong(pexp);

trL = size(p{1}{1},1);
ttS = chop([-preFrames+1:trL-preFrames]./frRateHz,2);

trLLong = size(pupilLong{1}{1},1);
ttLongS = chop([-preFrames+1:trLLong-preFrames]./frRateHz,2);

delayType = ms(1).audDelayType;

pszfig = figure;colormap gray;suptitle('pupil size')
hzfig = figure;colormap gray;suptitle('pupil horizontal pos')
vtfig = figure;colormap gray;suptitle('pupil vertical pos')
pupillongfig = figure;colormap gray
for idelay = 1:ntrtype
    pLong = [];
    hzLong = [];
    vtLong = [];
    for iexp = 1:sum(pexp)
        figure(pszfig)
        subplot(2,4,idelay)
        hold on
        plot(ttS,p{iexp}{idelay})
        title(delayType{idelay})
        figure(hzfig)
        subplot(2,4,idelay)
        hold on
        plot(ttS,hz{iexp}{idelay})
        title(delayType{idelay})
        figure(vtfig)
        subplot(2,4,idelay)
        hold on
        plot(ttS,vt{iexp}{idelay})
        title(delayType{idelay})
        
        if idelay <= ntrtype-2
            pLong = cat(2,pLong,pupilLong{iexp}{idelay});
            hzLong = cat(2,hzLong,horizLong{iexp}{idelay});
            vtLong = cat(2,vtLong,vertLong{iexp}{idelay});
        end
    end
    pLong = mean(pLong,2);
    hzLong = mean(hzLong,2);
    vtLong = mean(vtLong,2);
    
    figure(pszfig)
    figXAxis([],'time (s)',[])
    figYAxis([],'% change',[])
    figAxForm([])
    figure(hzfig)
    figXAxis([],'time (s)',[])
    figYAxis([],'diff',[])
    figAxForm([])
    figure(vtfig)
    figXAxis([],'time (s)',[])
    figYAxis([],'diff',[])
    figAxForm([])
    if idelay <= ntrtype-2
        figure(pupillongfig)
        hold on
        h = plot(ttLongS,pLong);
        h.Color = delay_colors(idelay,:);
        h.LineWidth = 2;
        delayS = str2num(ms(1).audDelayType{idelay})/1000;
        vline(delayS,'k--')
    end
end
figure(pupillongfig)
figXAxis([],'time(s)',[ttLongS(1) ttLongS(end)]);
figYAxis([],'% change',[])
figax = gca;
figax.TickDir = 'out';
figax.Box = 'off';
title('pupil size')
legend(ms(1).audDelayType(1:ntrtype-2))
print(fullfile(fn,'pupilsize_longTC'),'-dpdf','-fillpage');

figure(pszfig)
print(fullfile(fn,'pupilsize_AVdelay'),'-dpdf','-fillpage');
figure(hzfig)
print(fullfile(fn,'pupilhzpos_AVdelay'),'-dpdf','-fillpage');
figure(vtfig)
print(fullfile(fn,'pupilvtpos_AVdelay'),'-dpdf','-fillpage');



%% speed
sexp = false(1,nexp);
sexpt_str = cell(1,nexp);
sp = cell(1,nexp);
speedLong = cell(nexp,ntrtype-2);
for iexp = 1:nexp
    if ~isempty(ms(iexp).trSpeed)
        pw = ms(iexp).pre_event_wheel;
        combotrials = cellfun(@(x) cell2mat(x'),ms(iexp).trSpeed, 'unif', 0);
        bl = cellfun(@(x) mean(x(:,1:pw),2), combotrials, 'unif', 0);
        normsp = cellfun(@(x,y) bsxfun(@minus, x, y), combotrials, bl, 'unif', 0);
        sp{iexp} = cellfun(@(x) mean(x,1),normsp, 'unif', 0);
        
        sexp(iexp) = true;
        expt_str{iexp} = [ms(iexp).subnum '-' ms(iexp).date];
    end
    
    speeddata = ms(iexp).trSpeedLong;
    if ~isempty(speeddata)
        for idelay = 1:ntrtype-2
            speeddatadelay = speeddata{idelay};
            maxtriallength = max(cellfun(@length,speeddatadelay));        
            padtrialnan = cellfun(@(x) nan(1,maxtriallength-length(x)), speeddatadelay,'unif', 0);
            speeddatapadded = cellfun(@(x,y) cat(2,x,y), speeddatadelay, padtrialnan, 'unif', 0);
            combotrials = cell2mat(speeddatapadded');
            normsp = bsxfun(@minus, combotrials, mean(combotrials(:,1:pw),2));
%             speedLong{iexp,idelay} = mean(normsp,1);
            speedLong{iexp,idelay} = normsp;
        end
    end
end
expttriallengthmax = max(cellfun(@(x) size(x,2), speedLong(sexp,:)),[],1);
meanSpeedLong = cell(1,ntrtype-2);
errSpeedLong = cell(1,ntrtype-2);
for idelay = 1:ntrtype-2    
    if any(cellfun(@(x) size(x,2),speedLong(sexp,idelay)) ~= expttriallengthmax(idelay))
        padtrialnan = cellfun(@(x) nan(size(x,1),expttriallengthmax(idelay)-size(x,2)), speedLong(sexp,idelay),'unif', 0);
        speedLong(sexp,idelay) = cellfun(@(x,y) cat(2,x,y), speedLong(sexp,idelay), padtrialnan, 'unif', 0);
    end
    meanSpeedLong{idelay} = nanmean(cell2mat(speedLong(sexp,idelay)),1);
    errSpeedLong{idelay} = ste(cell2mat(speedLong(sexp,idelay)),1);
end
        
 
speed = sp(sexp);
preWheel = double(unique(cell2mat({ms(sexp).pre_event_wheel})));
wheelRateHz = double(unique(cell2mat({ms(sexp).wheelRateHz})));

trL = length(speed{1}{1});
ttS = chop([-preWheel+1:trL-preWheel]./wheelRateHz,2);
delayType = ms(1).audDelayType;

figure;colormap gray;suptitle('running speed')
for idelay = 1:ntrtype
    for iexp = 1:sum(sexp)
        subplot(2,4,idelay)
        hold on
        plot(ttS,speed{iexp}{idelay})
    end
    title(delayType{idelay})
    figXAxis([],'time (s)',[])
    figYAxis([],'speed (m/s)',[])
    figAxForm([])
end
print(fullfile(fn,'speed_AVdelay'),'-dpdf','-fillpage');

speed = cat(1,speed{1});
speed_all = cell(1,ntrtype);
for idelay = 1:ntrtype
    speed_all{idelay} = mean(cell2mat(speed(:,idelay)),1);
end

figure;
for idelay = 1:ntrtype
    h = plot(ttS,speed_all{idelay});
    h.Color = delay_colors(idelay,:);
    h.LineWidth = 1;
    hold on
end
figXAxis([],'time (s)',[])
figYAxis([],'speed (m/s)',[])
figAxForm([])
title(sprintf('mean speed across mice (n=%s)',num2str(sum(sexp))))
legend(cat(2,trType_str(1:ntrtype)),'location','northeastoutside')

speedlongfig = figure;
for idelay = 1:ntrtype-2
    trLLong = size(meanSpeedLong{idelay},2);
    ttLongS = chop([-preWheel+1:trLLong-preWheel]./wheelRateHz,2);
    figure(speedlongfig)
    hold on
    h = plot(ttLongS(1:end-2),meanSpeedLong{idelay}(1:end-2));
    h.Color = delay_colors(idelay,:);
    h.LineWidth = 2;
%     h.edge(1).Color = [1 1 1];
%     h.edge(2).Color = [1 1 1];
    delayS = str2num(ms(1).audDelayType{idelay})/1000;
    vline(delayS,'k--')
end
figure(speedlongfig)
figXAxis([],'time(s)',[ttLongS(1) ttLongS(end)]);
figYAxis([],'m/s',[])
figax = gca;
figax.TickDir = 'out';
figax.Box = 'off';
title('speed')
legend(ms(1).audDelayType(1:ntrtype-2))
print(fullfile(fn,'speed_longTC'),'-dpdf','-fillpage');

%% combine datasets
nCellsExpt = zeros(1,size(ms,2));
delayTCOri = cell(nori,ntrtype);
delayTCOri_err = cell(nori,ntrtype);
delayOri_ttest = cell(nori,ntrtype);
delayTCOri_boot = cell(nori,ntrtype);
avRespMod = cell(1,nori);
delayNsCorr_expt = cell(nori,ntrtype+1,size(ms,2));
audOnlyTC = [];
audOnly_ttest = [];
delayTCLong = cell(nori,ntrtype-2);
for iexp = 1:size(ms,2)
    for itype = 1:ntrtype
        dt = ms(iexp).delayTrials{itype};
        nCellsExpt(iexp) = size(dt,2);
        if itype == ntrtype
            audOnlyTC = cat(2,audOnlyTC,mean(dt,3));
            pre = squeeze(mean(dt(pre_win,:,:),1));
            resp = squeeze(mean(dt(resp_win,:,:),1));
            audOnly_ttest = cat(1,audOnly_ttest,ttest(resp,pre,'dim',2,'tail','right'));
        end            
        for iori = 1:nori
            ind = ms(iexp).trOri{itype} == oris(iori);
            delayTCOri{iori,itype} = cat(2,delayTCOri{iori,itype},mean(dt(:,:,ind),3));
            boot_temp = zeros(size(dt,1),size(dt,2),nboot);
            for iboot = 1:nboot
                ind2 = randsample(find(ind),sum(ind),1); 
                boot_temp(:,:,iboot) = mean(dt(:,:,ind2),3);
            end
            delayTCOri_boot{iori,itype} = cat(2,delayTCOri_boot{iori,itype},boot_temp);
            delayTCOri_err{iori,itype} = cat(2,delayTCOri_err{iori,itype},ste(dt(:,:,ind),3));
            pt = mean(dt(pre_win,:,ind),1);
            rt = mean(dt(resp_win,:,ind),1);
            delayOri_ttest{iori,itype} = cat(2,delayOri_ttest{iori,itype},ttest(rt,pt,'dim',3,'tail','right'));
            delayNsCorr_expt{iori,itype,iexp} = corrcoef(squeeze(rt)');
            
            if ismember(itype,1:ntrtype-2)
                dtl = ms(iexp).delayLong{itype};
                delayTCLong{iori,itype} = cat(2,delayTCLong{iori,itype},mean(dtl(:,:,ind),3));
            end
            
        % cells modulated by auditory
            if itype == 1;
                av = mean(dt(resp_win,:,ind),1);
            elseif itype == visOnly_ind;
                v = mean(dt(resp_win,:,ind),1);
                avRespMod{iori} = cat(2,avRespMod{iori},ttest2(v,av,'dim',3,'tail','both'));
                iti = mean(dt(1:preFrames,:,ind),1);
                delayNsCorr_expt{iori,ntrtype+1,iexp} = corrcoef(squeeze(iti)');
            end
        end
    end
end
cumCellsExpt = cumsum(nCellsExpt);

delayNsCorr = squeeze(mean(cellfun(@(x) mean(x(tril(x,1) > 0)), delayNsCorr_expt),1));
delayNsCorr_err = squeeze(ste(cellfun(@(x) mean(x(tril(x,1) > 0)), delayNsCorr_expt),1));

delayRespOri = cellfun(@(x) mean(x(resp_win,:),1) - mean(x(pre_win,:),1), delayTCOri,'unif',0);
delayRespOri_err = cellfun(@(x) mean(x(resp_win,:),1) - mean(x(pre_win,:),1), delayTCOri_err,'unif',0);
delayRespOri_boot = cellfun(@(x) mean(x(resp_win,:,:),1) - mean(x(pre_win,:,:),1), delayTCOri_boot,'unif',0);


%% noise correlation
figure;
for iexp = 1:size(ms,2)
    h = errorbar(1:ntrtype+1,delayNsCorr(:,iexp),delayNsCorr_err(:,iexp),'o');
    h.MarkerFaceColor = expt_colors(iexp,:);
    h.MarkerEdgeColor = [1 1 1];
    hold on
end
figXAxis([],'aud stim delay (ms)',[0 ntrtype+2],1:ntrtype+1,cat(2,trType_str,{'iti'}))
figYAxis([],'mean noise corr',[0 1])
figAxForm([])
leg = legend(expt_str,'location','northeastoutside');
title({'correlated activity within a trial type';'mean & ste across orientations within an experiment'})
title(leg,'expt')
%% cell orientation tuning
nc = size(delayTCOri{1},2);
delayTuning = cell(1,ntrtype);
delayTuning_boot = cell(1,ntrtype);
delayTuning_err = cell(1,ntrtype);
delayTuningTtest = cell(1,ntrtype);
for itype = 1:ntrtype
    delayTuning{itype} = cell2mat(delayRespOri(:,itype));
    delayTuning_boot{itype} = cell2mat(delayRespOri_boot(:,itype));
    delayTuning_err{itype} = cell2mat(delayRespOri(:,itype));
    delayTuningTtest{itype} = cell2mat(delayOri_ttest(:,itype));
end
[delayOriMax, delayOriMax_ind] = cellfun(@(x) max(x,[],1),delayTuning,'unif',0);

% for vis only data, find if max response is also significant response
visOnlyMaxInd = delayOriMax_ind{visOnly_ind};
visOnlyTtest = delayTuningTtest{visOnly_ind};
visOnlyMaxSN = zeros(1,nc);
for icell = 1:nc
    tt = visOnlyTtest(:,icell);
   visOnlyMaxSN(icell) = tt(visOnlyMaxInd(icell));
end

visOnlyMaxInd_sub1 = visOnlyMaxInd-1;
visOnlyMaxInd_sub1(visOnlyMaxInd_sub1 == 0) = nori;

%% von mises fit tuning curves
orisFit = [oris 180];
orisSmooth = 0:1:180;

delayTuningBeforeFit = cellfun(@(x) cat(1,x,x(1,:)), delayTuning, 'unif', 0);
delayTuningBeforeFit_boot = cellfun(@(x) cat(1,x,x(1,:,:)), delayTuning_boot, 'unif', 0);

% get 90% confidence intervals for vis only data tuning fits
visOnlyData = delayTuningBeforeFit{visOnly_ind}';
visOnlyData_resamp = permute(delayTuningBeforeFit_boot{visOnly_ind},[2,1,3]);
if nboot == 1000
    try
        load(fullfile(fn,['awData_audMod' ds '_tuningFitReliability']));
    catch
        theta_90 = [];
    end
    if length(theta_90) ~= size(visOnlyData,1)
        doBoot = logical(input('Fit reliability not calculated for all cells. Do fits for all cells (0/1)?'));
        if doBoot
            [~, theta_dist, theta_90] = vonmisesReliableFit(visOnlyData, visOnlyData_resamp,orisFit,nboot);
            save(fullfile(fn,['awData_audMod' ds '_tuningFitReliability']),'theta_90','theta_dist');
        else 
            cmlvCellsExpt = cumsum(nCellsExpt);
            n = length(theta_90);
            theta_90_temp = theta_90;
            theta_dist_temp = theta_dist;
            ind = cmlvCellsExpt(cmlvCellsExpt == n)+1;
            visOnlyData = visOnlyData(ind:end,:);
            visOnlyData_resamp = visOnlyData_resamp(ind:end,:,:);
            [~, theta_dist, theta_90] = vonmisesReliableFit(visOnlyData, visOnlyData_resamp,orisFit,nboot);
            
            theta_90 = cat(2,theta_90_temp,theta_90);
            theta_dist = cat(2,theta_dist_temp,theta_dist);
            save(fullfile(fn,['awData_audMod' ds '_tuningFitReliability']),'theta_90','theta_dist');            
        end
    end
else
    try
        load(fullfile(fn,['awData_audMod' ds '_tuningFitReliability']));
        if length(theta_90) ~= size(visOnlyData,1)
            [~, theta_dist, theta_90] = vonmisesReliableFit(visOnlyData, visOnlyData_resamp,orisFit,nboot);
        end
    catch
        [~, theta_dist, theta_90] = vonmisesReliableFit(visOnlyData, visOnlyData_resamp,orisFit,nboot);
    end
end

figure;
histogram(theta_90)

% get fits on raw data for each delay
delayRSquare = cell(1,ntrtype);
delayTuningFits = cell(1,ntrtype);
for idelay = 1:ntrtype
    data_temp = delayTuningBeforeFit{idelay};
    R_square = NaN(1,nc);
    y_fit = zeros(length(orisSmooth),nc);
    for icell = 1:nc
        [b, k, R, u, ~,R_square(icell)] = miaovonmisesfit_ori(deg2rad(orisFit),data_temp(:,icell));
        y_fit(:,icell) = b+R.*exp(k.*(cos(2.*(deg2rad(orisSmooth)-u))-1));
    end
    delayRSquare{idelay} = R_square;
    delayTuningFits{idelay} = y_fit;
end

[fitTuningPeak, fitTuningPeak_ind] = max(delayTuningFits{visOnly_ind},[],1);

[~, delayFitPeak_ind] = cellfun(@(x) max(x,[],1), delayTuningFits, 'unif',0);

goodRS = delayRSquare{visOnly_ind} > 0.9 & delayRSquare{visOnly_ind} > 0.9;

rlblFit = theta_90 < 30;
%%
audDelayInterneuronAnalysis
%% individual cells with significant difference bw vis and vis+aud conditions (at pref ori)
avModPref = zeros(1,nc);
for icell = 1:nc
    pref_ind = visOnlyMaxInd(icell);
    avModPref(icell) = avRespMod{pref_ind}(icell);
end

%% plot example cells fits with raw data
exCells = randsample(find(rlblFit),9);
exC = exCells(1);

figure;
for icell = 1:9
    exC = exCells(icell);
    subplot(3,3,icell)
    tc1 = delayTuningBeforeFit{visOnly_ind}(:,exC);
    tc2 = delayTuningFits{visOnly_ind}(:,exC);
    tc3 = delayTuningBeforeFit{1}(:,exC);
    tc4 = delayTuningFits{1}(:,exC);
    plot(orisFit,tc1);
    hold on
    plot(orisSmooth, tc2);
    plot(orisFit,tc3);
    plot(orisSmooth,tc4);
    figXAxis([],'orientation (deg)',[-1 181], orisFit,orisFit)
    figYAxis([],'dF/F',[])
    figAxForm([])
    title({['cell #' num2str(exC)];['theta_90=' num2str(theta_90(exC))];['Rsquare=' num2str(delayRSquare{5}(exC))]})
%     title({['cell #' num2str(exC)];['Rsquare=' num2str(delayRSquare{5}(exC))]})
    if icell == 9
        legend({'raw vis only', 'fit vis only', 'raw vis+aud', 'fit vis+aud'})
    end
end
%% are fit tuning peaks correlated?
%plot scatter
resp_lim = [0 180];
figure;
suptitle({'tuning peak for curves fit to each trial type'})
for idelay = 1:6
    subplot(2,3,idelay)
    scatter(fitTuningPeak_ind(rlblFit),delayFitPeak_ind{idelay}(rlblFit),50,'k.');
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{5},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
% print(fullfile(fn,'scatter_maxResp_eaDelay'),'-dpdf','-fillpage');

%% long data
prefTCLong_eachCell = getPrefTC(delayTCLong,visOnlyMaxInd);
orthTCLong_eachCell = getOrthTC(delayTCLong,visOnlyMaxInd);
prefTCLong = cellfun(@(x) mean(x(:,rlblFit),2),prefTCLong_eachCell,'unif',0);
orthTCLong = cellfun(@(x) mean(x(:,rlblFit),2),orthTCLong_eachCell,'unif',0);

timeVectorS = ((1:length(prefTCLong{1}))-preFrames)/frRateHz;
setFigParams4Print('landscape');figure;
for idelay = 1:length(prefTCLong)
    delayS = str2num(ms(1).audDelayType{idelay})/1000;
    tc = prefTCLong{idelay};
    h = plot(timeVectorS,tc);
    h.Color = delay_colors(idelay,:);
    h.LineWidth = 2;
    leg(idelay) = h;
    hold on
    vline(delayS,'k--')
end
figXAxis([],'time(s)',[timeVectorS(1) timeVectorS(end)]);
figYAxis([],'dF/F',[])
fig = gca;
fig.TickDir = 'out';
fig.Box = 'off';
title(sprintf('response to preferred orientation, n = %s tuned cells',num2str(sum(rlblFit))))
legend(leg,ms(1).audDelayType(1:ntrtype-2))

print(fullfile(fn,'meanPref_longTC'),'-dpdf','-fillpage');

figure;
for idelay = 1:length(prefTCLong)
    delayS = str2num(ms(1).audDelayType{idelay})/1000;
    tc = orthTCLong{idelay};
    h = plot(timeVectorS,tc);
    h.Color = delay_colors(idelay,:);
    h.LineWidth = 2;
    leg(idelay) = h;
    hold on
    vline(delayS,'k--')
end
figXAxis([],'time(s)',[timeVectorS(1) timeVectorS(end)]);
figYAxis([],'dF/F',[])
fig = gca;
fig.TickDir = 'out';
fig.Box = 'off';
title(sprintf('response to orthogonal orientation, n = %s tuned cells',num2str(sum(rlblFit))))
legend(leg,ms(1).audDelayType(1:ntrtype-2))

print(fullfile(fn,'meanOrth_longTC'),'-dpdf','-fillpage');
%% get tuning curves aligned to max response

%shift curves to align max (to vis only) on position 4
fit_center_ind = find(orisSmooth == 90)-1;
center_ind = find(oris == 90)-1;
% oris_shift = circshift(orisSmooth,center_ind-1,2);
% oris_shift(1:center_ind-1) = -fliplr(orisSmooth(2:center_ind));
delayTuningFitsShift = cell(1,ntrtype);
delayTuningShift = cell(1,ntrtype);
for idelay = 1:ntrtype
    fit_tc_all = delayTuningFits{idelay};
    fit_tc_shift = zeros(size(fit_tc_all));
    tc_all = delayTuning{idelay};
    tc_shift = zeros(size(tc_all));
    for icell = 1:nc
        pref_ind = fitTuningPeak_ind(icell);
        fit_tc = fit_tc_all(:,icell);
        nshift = fit_center_ind-pref_ind;
        fit_tc_shift(:,icell) = circshift(fit_tc,nshift);
        
        ori_ind = round(fitTuningPeak_ind(icell)/30)*30;
        if ori_ind == 180
            ori_ind = 0;
        end
        pref_ind = find(oris == ori_ind);
        tc = tc_all(:,icell);
        nshift = center_ind-pref_ind;
        tc_shift(:,icell) = circshift(tc,nshift);
    end
    delayTuningFitsShift{idelay} = fit_tc_shift;
    delayTuningShift{idelay} = tc_shift;
end

% normalize raw data to vis only condition

delayTuningNorm = cellfun(@(x) bsxfun(@rdivide,x,delayOriMax{visOnly_ind}), delayTuningShift ,'unif',0);
delayTuningNorm = cellfun(@(x) cat(1,x(end,:),x), delayTuningNorm, 'unif', 0);
tuning_mean = cellfun(@(x) mean(x(:,logical(rlblFit)),2), delayTuningNorm, 'unif',0);
tuning_ste = cellfun(@(x) ste(x(:,logical(rlblFit)),2), delayTuningNorm, 'unif',0);

tuning_mean_expt = cell(size(ms,2),ntrtype);
tuning_ste_expt = cell(size(ms,2),ntrtype);
for iexp = 1:size(ms,2)
    if iexp == 1
        offset = 0;
    else
        offset = cumCellsExpt(iexp-1);
    end
    ind = intersect(find(rlblFit),1+offset:nCellsExpt(iexp)+offset);
    tuning_mean_expt(iexp,:) = cellfun(@(x) mean(x(:,ind),2), delayTuningNorm, 'unif',0);
    tuning_ste_expt(iexp,:) = cellfun(@(x) ste(x(:,ind),2), delayTuningNorm, 'unif',0);
end

% normalize fit curves to vis only condition
delayTuningFitsShiftNorm = cellfun(@(x) bsxfun(@rdivide,x,fitTuningPeak), delayTuningFitsShift ,'unif',0);
fitTuning_mean = cellfun(@(x) nanmean(x(:,logical(rlblFit)),2), delayTuningFitsShiftNorm, 'unif',0);
fitTuning_ste = cellfun(@(x) ste(x(:,logical(rlblFit)),2), delayTuningFitsShiftNorm, 'unif',0);
%% plot normalized time-course across tuned cells
% timing info
trL = size(delayTCOri{1,1},1);
ttS = chop([-preFrames+1:trL-preFrames]./frRateHz,2);
x_axis_label = [-1:0.5:1];
audStimTime = [-1:0.25:0];

%get tc for pref ori only and orthogonal, normalized to vis only mean
%response
fit_center_ind = find(orisSmooth == 90)-1;
tuningPeak = round(fitTuningPeak_ind/30)*30;
tuningPeak(tuningPeak == 180) = 0;
delayTCPrefNorm = cell(1,ntrtype);
delayTCPrefNorm_err = cell(1,ntrtype);
delayTCOrthNorm = cell(1,ntrtype);
delayTCOrthNorm_err = cell(1,ntrtype);
for idelay = 1:ntrtype
    tc = nan(trL,nc);
    for icell = 1:nc
        pref_ind = oris == tuningPeak(icell);
        orth_ind = circshift(pref_ind,4,2);
        tc = delayTCOri{pref_ind,idelay}(:,icell);
        visonlytc = delayTCOri{pref_ind,visOnly_ind}(:,icell);
        tc_err(:,icell) = delayTCOri_err{pref_ind,idelay}(:,icell);
        tco = delayTCOri{orth_ind,idelay}(:,icell);
        tco_err(:,icell) = delayTCOri{orth_ind,idelay}(:,icell);
        m = mean(visonlytc(resp_win,:),1);
%         me = max(tc_err,[],1);
        tc_norm(:,icell) = tc./m;
        tco_norm(:,icell) = tco./m;
%         tc_err_norm(:,icell) = tc_err./me;
    end
    delayTCPrefNorm{idelay} = tc_norm;
    delayTCPrefNorm_err{idelay} = tc_err;
    delayTCOrthNorm{idelay} = tco_norm;
    delayTCOrthNorm_err{idelay} = tco_err;
end

%get tc for pref ori only and orthogonal, all max response normalized to 1
delayTCPrefEq = cellfun(@(x) mean(x(:,rlblFit),2)./mean(mean(x(resp_win,rlblFit),2)), delayTCPrefNorm, 'unif',0);
delayTCOrthEq = cellfun(@(x) mean(x(:,rlblFit),2)./mean(mean(x(resp_win,rlblFit),2)), delayTCOrthNorm, 'unif',0);
%% example cell time-courses
exCells = randsample(find(rlblFit & avModPref),3);
figure;
for icell = 1:3
    exC = exCells(icell);
    iplot = 1+((icell-1)*2);
    subplot (3,2,iplot)
    for idelay = 1:ntrtype-1
        tc = mean(delayTCPrefNorm{idelay}(:,exC),2);
        tc_err = ste(delayTCPrefNorm{idelay}(:,exC),2);
        h = shadedErrorBar(ttS,tc,tc_err,'k');
        if idelay < visOnly_ind
            h.mainLine.Color = delay_colors(idelay,:);
        end
        h.mainLine.LineWidth = 1;
        hold on
        leg(idelay) = h.mainLine;
    end
    leg(idelay+1) = vline(0,'k--');
    figXAxis([],'time (s)',[-1.1 1.1],x_axis_label,x_axis_label)
    figYAxis([],'dF/F norm',[])
    figAxForm([])
    vline(audStimTime,'k--')
    title({'preferred orientation';['cell #' num2str(exC)]})
%     legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')
    iplot = icell*2;
    subplot (3,2,iplot)
    for idelay = 1:ntrtype-1
        tc = mean(delayTCOrthNorm{idelay}(:,exC),2);
        tc_err = ste(delayTCOrthNorm{idelay}(:,exC),2);
        h = shadedErrorBar(ttS,tc,tc_err,'k');
        if idelay < visOnly_ind
            h.mainLine.Color = delay_colors(idelay,:);
        end
        h.mainLine.LineWidth = 1;
        hold on
        leg(idelay) = h.mainLine;
    end
    leg(idelay+1) = vline(0,'k--');
    figXAxis([],'time (s)',[-1.1 1.1],x_axis_label,x_axis_label)
    figYAxis([],'dF/F norm',[])
    figAxForm([])
    vline(audStimTime,'k--')
    title({'orthogonal orientation';['cell #' num2str(exC)]})
%     legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')
end

print(fullfile(fn,'exCells_normTCprefstim_AVdelay'),'-dpdf','-fillpage');

%% all tuned cells time-course, normalized to vis-only
figure;
subplot 121
for idelay = 1:ntrtype-1
    tc = mean(delayTCPrefNorm{idelay}(:,rlblFit),2);
    tc_err = ste(delayTCPrefNorm{idelay}(:,rlblFit),2);
    h = shadedErrorBar(ttS,tc,tc_err,'k');
    if idelay < visOnly_ind
        h.mainLine.Color = delay_colors(idelay,:);
    end
    h.mainLine.LineWidth = 1;
    hold on
    leg(idelay) = h.mainLine;
end
leg(idelay+1) = vline(0,'k--');
figXAxis([],'time from vis stim onset (1s stim) (s)',[-1.1 1.1],x_axis_label,x_axis_label)
figYAxis([],'dF/F normalized to vis only resp',[-0.2 2])
figAxForm([])
vline(audStimTime,'k--')
title('response to preferred orientation (tuned neurons)')
legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')
subplot 122
for idelay = 1:ntrtype-1
    tc = mean(delayTCOrthNorm{idelay}(:,rlblFit),2);
    tc_err = ste(delayTCOrthNorm{idelay}(:,rlblFit),2);
    h = shadedErrorBar(ttS,tc,tc_err,'k');
    if idelay < visOnly_ind
        h.mainLine.Color = delay_colors(idelay,:);
    end
    h.mainLine.LineWidth = 1;
    hold on
    leg(idelay) = h.mainLine;
end
leg(idelay+1) = vline(0,'k--');
figXAxis([],'time from vis stim onset (1s stim) (s)',[-1.1 1.1],x_axis_label,x_axis_label)
figYAxis([],'dF/F normalized to vis only resp',[-0.2 2])
figAxForm([])
vline(audStimTime,'k--')
title('response to orthogonal orientation (tuned neurons)')
legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')

print(fullfile(fn,'norm2visonly_TCprefstim_AVdelay'),'-dpdf','-fillpage');
%% all tuned cells time-course, normalized to 1
figure;
subplot 121
for idelay = 1:ntrtype-1
    tc = delayTCPrefEq{idelay};
%     tc_err = ste(delayTCPrefEq{idelay}(:,rlblFit),2);
    h = plot(ttS,tc,'k');
    if idelay < visOnly_ind
        h.Color = delay_colors(idelay,:);
    end
    h.LineWidth = 1;
    hold on
    leg(idelay) = h;
end
leg(idelay+1) = vline(0,'k--');
figXAxis([],'time from vis stim onset (1s stim) (s)',[-0.5 1],x_axis_label,x_axis_label)
figYAxis([],'dF/F norm to mean resp',[-.3 1.5])
figAxForm([])
vline(audStimTime,'k--')
title('response to preferred orientation (tuned neurons)')
legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')
subplot 122
for idelay = 1:ntrtype-1
    tc = delayTCOrthEq{idelay};
    h = plot(ttS,tc,'k');
    if idelay < visOnly_ind
        h.Color = delay_colors(idelay,:);
    end
    h.LineWidth = 1;
    hold on
    leg(idelay) = h;
end
leg(idelay+1) = vline(0,'k--');
figXAxis([],'time from vis stim onset (1s stim) (s)',[-0.5 1],x_axis_label,x_axis_label)
figYAxis([],'dF/F norm to mean resp',[-.3 1.5])
figAxForm([])
vline(audStimTime,'k--')
title('response to orthogonal orientation (tuned neurons)')
legend(leg,cat(2,trType_str(1:ntrtype-1),'aud stim delay'),'location','northwest')

print(fullfile(fn,'normTCprefstim_AVdelay'),'-dpdf','-fillpage');

%% plot normalized curves across tuned cells

figure;
subplot 131
plot(orisSmooth,delayTuningFitsShiftNorm{visOnly_ind}(:,logical(rlblFit)))
figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
figYAxis([],'norm dF/F',[-1 2])
title('vis only')

subplot 132
plot(orisSmooth,delayTuningFitsShiftNorm{1}(:,logical(rlblFit)))
figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
figYAxis([],'norm dF/F',[-1 2])
title('vis+aud')

subplot 133
h = errorbar(orisSmooth,fitTuning_mean{visOnly_ind},fitTuning_ste{visOnly_ind},'k');
h.LineWidth = 1;
hold on
h = errorbar(orisSmooth,fitTuning_mean{1},fitTuning_ste{1},'m');
h.LineWidth = 1;
figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
figYAxis([],'norm dF/F',[0 1.2])
title('mean across responsive cells')
legend({'vis only';'vis+aud'})

print(fullfile(fn,'normTuning_AV'),'-dpdf','-fillpage');

figure;
for idelay = 1:ntrtype-2
    subplot(1,ntrtype-2,idelay)
    h = errorbar(orisSmooth,fitTuning_mean{visOnly_ind},fitTuning_ste{visOnly_ind},'k');
    h.LineWidth = 1;
    hold on
    h = errorbar(orisSmooth,fitTuning_mean{idelay},fitTuning_ste{idelay},'m');
    h.LineWidth = 1;
    figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
    figYAxis([],'norm dF/F',[0 1.2])
    title([trType_str{idelay} 'ms delay'])
    legend({'vis only';'vis+aud'})
end

print(fullfile(fn,'normTuningFits_AVdelay'),'-dpdf','-fillpage');

%% raw data aligned tuning
tuning_mean_fit = cell(1,ntrtype);
for idelay = 1:ntrtype
    [b, k, R, u, ~,R_square(icell)] = miaovonmisesfit_ori(deg2rad(orisFit),tuning_mean{idelay});
    y_fit = b+R.*exp(k.*(cos(2.*(deg2rad(orisSmooth)-u))-1));
    tuning_mean_fit{idelay} = y_fit;
end
setFigParams4Print('landscape');figure;
for idelay = 1:ntrtype-2
    subplot(1,ntrtype-2,idelay)
    h = errorbar(orisFit,tuning_mean{visOnly_ind},tuning_ste{visOnly_ind},'k');
    h.LineWidth = 2;
    hold on
    h = errorbar(orisFit,tuning_mean{idelay},tuning_ste{idelay},'m');
    h.LineWidth = 2;
    figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
    figYAxis([],'norm dF/F',[0 1.3])
    title([trType_str{idelay} 'ms delay'])
    if idelay == 1
        legend({'vis only';'vis+aud'},'location','northwest')
    end
    figAxForm([])
    
%     subplot(2,4,idelay+4)
%     h = plot(orisSmooth,tuning_mean_fit{5},'k');
%     h.LineWidth = 2;
%     hold on
%     h = plot(orisSmooth,tuning_mean_fit{idelay},'m');
%     h.LineWidth = 2;
%     figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
%     figYAxis([],'fit dF/F',[0 1.3])
%     title([trType_str{idelay} 'ms delay'])
%     figAxForm([])
end

print(fullfile(fn,'normTuning_AVdelay'),'-dpdf','-fillpage');

%individual experiments
for iexp = 1:size(ms,2)
    setFigParams4Print('landscape');figure;
    if iexp == 1;
        offset = 0;
        ncells = sum(rlblFit(1:cumCellsExpt(iexp)));
    else
        offset = cumCellsExpt(iexp-1);
        ncells = sum(rlblFit(1+offset:cumCellsExpt(iexp)));
    end
    expt_name = [ms(iexp).subnum '-' ms(iexp).date];
    suptitle({expt_name;[num2str(ncells) ' cells']})
    for idelay = 1:ntrtype-2
        subplot(1,ntrtype-2,idelay)
        h = errorbar(orisFit,tuning_mean_expt{iexp,visOnly_ind},tuning_ste_expt{visOnly_ind},'k');
        h.LineWidth = 2;
        hold on
        h = errorbar(orisFit,tuning_mean_expt{iexp,idelay},tuning_ste_expt{idelay},'m');
        h.LineWidth = 2;
        figXAxis([],'orientation (deg)',[-1 181],orisFit,orisFit)
        figYAxis([],'norm dF/F',[0 1.3])
        title([trType_str{idelay} 'ms delay'])
        if idelay == 1
            legend({'vis only';'vis+aud'},'location','northwest')
        end
        figAxForm([])
    end
    
    print(fullfile(fn,'expt sum',[expt_name '_normTuning_AVdelay']),'-dpdf','-fillpage');
end

%% plot example cells orientation curves for each delay
exCells = randsample(find(rlblFit),16);
exC = exCells(1);


figure;
suptitle(['cell #' num2str(exC)])
subplot 121
for itype = 1:ntrtype
    tc = delayTuning{itype}(:,exC);
    err = delayTuning_err{itype}(:,exC);
    hold on
    h = errorbar(oris,tc,err,'ko-');
    h.Color = delay_colors(itype,:);
    h.LineWidth = 1;
    h.MarkerFaceColor = delay_colors(itype,:);
end
legend(trType_str,'location','northeastoutside')
figXAxis([],'orientation (deg)',[oris(1)-10 oris(end)+10],oris,oris)
figYAxis([],'dF/F',[])
figAxForm([])
hold on
plot(oris(visOnlyMaxInd(exC)),tc(visOnlyMaxInd(exC))+0.1,'k*')
title('average responses')

subplot 122
for itype = 1:ntrtype
    tc = delayTuningFits{itype}(:,exC);
    hold on
    h = plot(orisSmooth,tc,'ko-');
    h.Color = delay_colors(itype,:);
    h.LineWidth = .5;
    h.MarkerFaceColor = delay_colors(itype,:);
end
legend(trType_str,'location','northeastoutside')
figXAxis([],'orientation (deg)',[oris(1)-10 oris(end)+10],oris,oris)
figYAxis([],'dF/F',[])
figAxForm([])
hold on
plot(oris(visOnlyMaxInd(exC)),tc(visOnlyMaxInd(exC))+0.1,'k*')
title('von mises fits')

print(fullfile(fn,'exCell_tuning'),'-dpdf','-fillpage');

%% plot example cell time-course for each delay, multiple orientations
trLenFr = size(delayTCOri{1,1},1);
trMs = round((-preFrames+1:trLenFr-preFrames)/frRateHz*1000);
exCells = randsample(find(rlblFit),16);
exC = exCells(1);
figure;
suptitle(['example cell #' num2str(exC)])
for iori = 1:nori
    tc1 = delayTCOri{iori,1}(:,exC);
    err1 = delayTCOri_err{iori,1}(:,exC);
    tc2 = delayTCOri{iori,visOnly_ind}(:,exC);
    err2 = delayTCOri_err{iori,visOnly_ind}(:,exC);
    subplot(1,nori,iori)
    h = shadedErrorBar(trMs,tc1,err1,'m-');
    leg(1) = h.mainLine;
    hold on
    h = shadedErrorBar(trMs,tc2,err2,'k-');
    leg(2) = h.mainLine;
    hold on
    vline(0,'k--')
    figXAxis([],'time from vis stim on (ms)',[-250 1000])
    figYAxis([],'dF/F',[-0.05 0.15])
    figAxForm([])
    title([num2str(oris(iori)) ' deg'])
    if iori == 1
        legend(leg,{'vis+aud';'vis only'},'location','northwest');
    end
end

print(fullfile(fn,'exCell_tc_tuning'),'-dpdf','-fillpage');

% % avTuningB4Fit = cat(1,delayTuning{1},delayTuning{1}(1,:))';
% % vTuningB4Fit = cat(1,delayTuning{visOnly_ind},delayTuning{visOnly_ind}(1,:))';
% % avTuningB4Fit_boot = permute(cat(1,delayTuning_boot{1},delayTuning_boot{1}(1,:,:)),[2,1,3]);
% % vTuningB4Fit_boot = permute(cat(1,delayTuning_boot{visOnly_ind},delayTuning_boot{visOnly_ind}(1,:,:)),[2,1,3]);
% % avFits = vonmisesReliableFit(avTuningB4Fit, avTuningB4Fit_boot,orisFit,nboot);
% % vFits = vonmisesReliableFit(vTuningB4Fit, vTuningB4Fit_boot,orisFit,nboot);
% % 
% % avFitMod_h = nan(1,nc);
% % for icell = 1:nc
% %     pref_ind = fitTuningPeak_ind(icell);
% %     avFitMod_h(icell) = ttest(squeeze(vFits(pref_ind,:,icell)),squeeze(avFits(pref_ind,:,icell)),'tail','both');
% % end
% % 
% % avRespMod_h = nan(1,nc);
% % for icell = 1:nc
% %     pref_ind = visOnlyMaxInd(icell);
% %     avRespMod_h(icell) = avRespMod{pref_ind}(icell);
% % end
%% plot scatter of responses to each delay type compared to vis only 
% (prefferred orienation)
fitTuningOrth = fitTuningPeak_ind-90;
fitTuningOrth(fitTuningOrth<0) = fitTuningPeak_ind(fitTuningOrth<0)+90;
fitTuningOrth = fitTuningOrth+1;
respTuningOrth = visOnlyMaxInd-3;
respTuningOrth(respTuningOrth<1) = visOnlyMaxInd(respTuningOrth<1)+3;
% get response in pref ori, only for signif cells
visOnlyFit_pref = NaN(1,nc);
visOnlyFit_orth = NaN(1,nc);
delayFit_pref = cell(1,ntrtype-2);
delayFit_orth = cell(1,ntrtype-2);
visOnlyResp_pref = NaN(1,nc);
visOnlyResp_orth = NaN(1,nc);
delayResp_pref = cell(1,ntrtype-2);
delayResp_orth = cell(1,ntrtype-2);
for icell = 1:nc
    %fit peak
    ind = fitTuningPeak_ind(icell);
    ind2 = fitTuningOrth(icell);

    visOnlyFit_pref(icell) = delayTuningFits{visOnly_ind}(ind,icell);
    visOnlyFit_orth(icell) = delayTuningFits{visOnly_ind}(ind2,icell);
    for idelay = 1:ntrtype-2
        delayFit_pref{idelay}(icell) = delayTuningFits{idelay}(ind,icell);
        delayFit_orth{idelay}(icell) = delayTuningFits{idelay}(ind2,icell);
    end
    % raw data peak
    ind = visOnlyMaxInd(icell);
    ind2 = respTuningOrth(icell);

    visOnlyResp_pref(icell) = delayTuning{visOnly_ind}(ind,icell);
    visOnlyResp_orth(icell) = delayTuning{visOnly_ind}(ind2,icell);
    for idelay = 1:ntrtype-2
        delayResp_pref{idelay}(icell) = delayTuning{idelay}(ind,icell);
        delayResp_orth{idelay}(icell) = delayTuning{idelay}(ind2,icell);
    end
end

%plot scatter expt identification
expt_ind = ones(1,nc);
for iexp = 1:size(ms,2)
    if iexp == 1
        offset = 0;
    else
        offset = offset+nCellsExpt(iexp-1);
        expt_ind(1+offset:nCellsExpt(iexp)+offset) = iexp;
    end
end
%%  raw data scatter  
resp_lim = [-0.1 1];
figure;
suptitle('response to preferred orientation (raw data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
%     h = scatter(visOnlyResp_pref(rlblFit),delayResp_pref{idelay}(rlblFit),20,'o');
%     h.MarkerEdgeColor = [0.5 0.5 0.5];
    for iexp = 1:size(ms,2)
        hold on
        h = scatter(visOnlyResp_pref(rlblFit & avModPref & expt_ind == iexp),delayResp_pref{idelay}(rlblFit & avModPref & expt_ind == iexp),50,'o');
        h.MarkerFaceAlpha = .75;
        h.MarkerFaceColor = expt_colors(iexp,:);
        h.MarkerEdgeColor = [1 1 1];
    end
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{visOnly_ind},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
print(fullfile(fn,'scatter_prefResp_eaDelay'),'-dpdf','-fillpage');
resp_lim = [-0.01 0.11];
figure;
suptitle('response to orthogonal orientation (raw data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
%     h = scatter(visOnlyResp_pref(rlblFit),delayResp_pref{idelay}(rlblFit),20,'o');
%     h.MarkerEdgeColor = [0.5 0.5 0.5];
    for iexp = 1:size(ms,2)
        hold on
        h = scatter(visOnlyResp_orth(rlblFit & avModPref & expt_ind == iexp),delayResp_orth{idelay}(rlblFit & avModPref & expt_ind == iexp),50,'o');
        h.MarkerFaceAlpha = .75;
        h.MarkerFaceColor = expt_colors(iexp,:);
        h.MarkerEdgeColor = [1 1 1];
    end
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{visOnly_ind},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
legend(expt_str,'location','southeast')
print(fullfile(fn,'scatter_orthResp_eaDelay'),'-dpdf','-fillpage');

%%  fit data scatter  
resp_lim = [-0.05 0.3];
figure;
suptitle('response to preferred orientation (fit data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
%     h = scatter(visOnlyResp_pref(rlblFit),delayResp_pref{idelay}(rlblFit),20,'o');
%     h.MarkerEdgeColor = [0.5 0.5 0.5];
    for iexp = 1:size(ms,2)
        hold on
        h = scatter(visOnlyFit_pref(rlblFit & expt_ind == iexp),delayFit_pref{idelay}(rlblFit & expt_ind == iexp),50,'o');
        h.MarkerFaceAlpha = .75;
        h.MarkerFaceColor = expt_colors(iexp,:);
        h.MarkerEdgeColor = [1 1 1];
    end
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{visOnly_ind},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
print(fullfile(fn,'scatter_prefFit_eaDelay'),'-dpdf','-fillpage');
resp_lim = [-0.01 0.11];
figure;
suptitle('response to orthogonal orientation (fit data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
%     h = scatter(visOnlyResp_pref(rlblFit),delayResp_pref{idelay}(rlblFit),20,'o');
%     h.MarkerEdgeColor = [0.5 0.5 0.5];
    for iexp = 1:size(ms,2)
        hold on
        h = scatter(visOnlyFit_orth(rlblFit & expt_ind == iexp),delayFit_orth{idelay}(rlblFit & expt_ind == iexp),50,'o');
        h.MarkerFaceAlpha = .75;
        h.MarkerFaceColor = expt_colors(iexp,:);
        h.MarkerEdgeColor = [1 1 1];
    end
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{visOnly_ind},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
legend(expt_str,'location','southeast')
print(fullfile(fn,'scatter_orthFit_eaDelay'),'-dpdf','-fillpage');

%% distribution of responses
resp_lim = [-0.02 0.15];
figure;
suptitle('response to preferred orientation (fit data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
    h = cdfplot(visOnlyFit_pref(rlblFit));
    h.Color = 'k';
    h.LineWidth = 1;
    hold on
    h = cdfplot(delayFit_pref{idelay}(rlblFit));  
    h.Color = 'm';
    h.LineWidth = 1;
    figXAxis([],'dF/F',resp_lim)
    figYAxis([],'% cells',[0 1])
    figAxForm([])
    title([trType_str{idelay} 'ms delay'])
end
legend({'vis only';'vis+aud'},'location','southeast')
print(fullfile(fn,'cdf_prefFit_eaDelay'),'-dpdf','-fillpage');
figure;
suptitle('response to orthogonal orientation (fit data)')
for idelay = 1:ntrtype-2
    subplot(2,3,idelay)
    h = cdfplot(visOnlyFit_orth(rlblFit));
    h.Color = 'k';
    h.LineWidth = 1;
    hold on
    h = cdfplot(delayFit_orth{idelay}(rlblFit));  
    h.Color = 'm';
    h.LineWidth = 1;
    figXAxis([],'dF/F',resp_lim)
    figYAxis([],'% cells',[0 1])
    figAxForm([])
    title([trType_str{idelay} 'ms delay'])
end
legend({'vis only';'vis+aud'},'location','southeast')
print(fullfile(fn,'cdf_orthFit_eaDelay'),'-dpdf','-fillpage');

%% running correlation
speedRespFig = figure;
suptitle('response to preferred (raw data) by running speed')
speed_ind = [];
delaySpeedTC = cell(1,ntrtype);
delayDsTC = cell(1,ntrtype);
delayDsTCErr = cell(1,ntrtype);
delayPreSpeed = cell(1,ntrtype);
delayRespSpeed = cell(1,ntrtype);
delayRunTtest = cell(1,ntrtype);
for iexp = 1:size(ms,2)
    for idelay = 1:ntrtype
        dt = ms(iexp).delayTrials{idelay};
        if isempty(ms(iexp).trSpeed)
            if idelay == 1
                speed_ind = cat(2,speed_ind,zeros(1,size(dt,2)));
            end
        else
            preWheel = ms(iexp).pre_event_wheel;
            whRateHz = ms(iexp).wheelRateHz;
            sp_resp_win = preWheel + round(.1*whRateHz):preWheel + round((respWinS+.1)*whRateHz);
            sp_pre_win = preWheel-round(respWinS*whRateHz)+1:preWheel;
            sp = ms(iexp).trSpeed{idelay};
            s = mean(sp(:,sp_pre_win),2);
            ds = bsxfun(@minus, sp,s);
            rsp =  squeeze(mean(dt(resp_win,:,:),1));
            ntr = size(rsp,2);
            pref_ind = visOnlyMaxInd(expt_ind == iexp);
            rsp_pref = nan(1,ntr);
            for iori = 1:nori
                tr_ind = ms(iexp).trOri{idelay} == oris(iori);
                cell_ind = pref_ind == iori;
               rsp_pref(tr_ind) = mean(rsp(cell_ind,tr_ind),1);
            end
            delaySpeedTC{idelay} = cat(1,delaySpeedTC{idelay},mean(sp,1));
            delayDsTC{idelay} = cat(1,delayDsTC{idelay},mean(ds,1));
            delayDsTCErr{idelay} = cat(1,delayDsTCErr{idelay},ste(ds,1));
            delayPreSpeed{idelay} = cat(1,delayPreSpeed{idelay},s);
            delayRespSpeed{idelay} = cat(1,delayRespSpeed{idelay},mean(sp(:,sp_resp_win),2));
            run_ttest = ttest(mean(sp(:,sp_resp_win),2),s);
            delayRunTtest{idelay} = cat(1,delayRunTtest{idelay},run_ttest);
            if idelay == 1
                speed_ind = cat(2,speed_ind,ones(1,size(dt,2)));
            end
%             figure(speedRespFig)
%             subplot(4,2,idelay)
%             hold on
%             h = scatter(mean(sp(:,sp_resp_win),2),rsp_pref,30,'o');
%             h.MarkerFaceAlpha = .75;
%             h.MarkerFaceColor = expt_colors(iexp,:);
%             h.MarkerEdgeColor = [1 1 1];
        end 
%         if iexp == size(ms,2)
%             x_lim = [-2 5];
%             y_lim = [-0.05 0.1];
%             title(trType_str{idelay})
%             figXAxis([],'speed (m/s)',[])
%             figYAxis([],'dF/F',[])
%             figAxForm([])
%         end
    end
end
figure(speedRespFig)
print(fullfile(fn,'runResp2Pref'),'-dpdf','-fillpage');



ttP = double(-preWheel+1:length(delayDsTC{1})-preWheel);
ttS = ttP/double(whRateHz);

figure;
leg = [];
for idelay = 1:ntrtype
    hold on
    h = shadedErrorBar(ttS,delayDsTC{idelay},delayDsTCErr{idelay},[],1);
    h.mainLine.Color = delay_colors(idelay,:);
    h.mainLine.LineWidth = 2;
    leg(idelay) = h.mainLine;
end
figXAxis([],'time from trial start (s)',[])
figYAxis([],'normalized speed (m/s)',[])
figAxForm([])
legend(leg,trType_str,'location','southeastoutside')
title('speed normalized to pre-response window by subtraction')
print(fullfile(fn,'runTC'),'-dpdf','-fillpage');

sp_lim = [-2 2];
figure;
for idelay = 1:ntrtype
    hold on
    sp = delayRespSpeed{idelay} - delayPreSpeed{idelay};
    h = cdfplot(sp);
    h.Color = delay_colors(idelay,:);
    h.LineWidth = 2;
end
figXAxis([],'speed (m/s)',sp_lim)
figYAxis([],'% trials',[0 1])
figAxForm([])
legend(trType_str,'location','southeastoutside')
title('running speed (resp win - base win)')


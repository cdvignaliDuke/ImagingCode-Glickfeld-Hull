%% full size-curve fitting + bootstrap analysis code
% adapted from sizeCurveFittingSmooth + bootstrapTest
% updated 4/23/18 KM

% loads size tuning data from size x contrast experiments along with retinotopy fit results
% will load existing or create new sizeTuneData.m (dF/F responses + cellDists)
% then runs fitting with both of the single- or double-sigmoid models:
% Model1:
% Model2:
% using fminsearchbnd with a modified objective function including a smoothness penalty
% subsequently runs bootstrapping on highest contrast condition
% will skip and load results if already exist for fitting and bootstrapping
% Results:
% presents individual cells with fits at each condition and prefSize+SI
% presents overall prefSize, SI, RF-stim plots

%% Load data
% currently need to run sizeTuningAfterRet, runs through only one dataset
%
% here load/define variables from chosen data set including:
% time course data
% size conditions, indices
% retinotopy fits, specifically RF centers
% good fit indices

% then manually define stim Az+El
% select only cells with RF center within range of stim (<7 deg)
% define size tuning curve
% ISSUE: load contrasts, need to have saved aside from loading input
% size tuning data - sizeTune and sizeSEM as (size, con, cell)

clear all;clc;
%i840 %180505 %2002
mouse = 'i884';
date = '180718';
ImgFolder = char('003');
time = char('2035');

RetImgFolder = char('002');
override_all = 0; % override_all to set override for each section (0=load previous, 1=write over)

nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

fprintf(['Size tuning curve analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

% load behavior/experimental data input
fprintf(['Loading input from size tuning runs: ' run_str '\n'])
fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(1,:) '.mat'];
load(fName);
for i=1:length(input.tGratingContrast) % replace int64(con=1) with double
    if ~(class(input.tGratingContrast{i})=="double")
        input.tGratingContrast{i} = double(input.tGratingContrast{i});
    end
end

% load tc data
% loads 'tc_dfof', 'tuning_mat', 'szs', 'Ind_struct'
% tc_dfof is timecourse data - (trialFrames x nCells x nTrials)
% szs is list of size conditions
% Ind_struct is 1xlength(szs) struct with trial indices for each size condition
fprintf(['Loading timecourses from size tuning runs: ' run_str '\n'])
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']))

% load experimental conditions from input
% trial variables
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2);
% contrast trial indices and range
conTrials = cell2mat(input.tGratingContrast);
cons = unique(conTrials);
nCon = length(cons);
% make conInds, a cell array with vectors of indices for each contrast
conInds = cell(nCon,1);
for i = 1:nCon
    conInds{i} = find(conTrials == cons(i));
end
nSize = length(szs);

% load retinotopy data
nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);

% loads 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind'
% lbub_fits contains RF fit data, extract RF centers from here
% goodfit_ind
fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
load(fn_out);
cellAz = lbub_fits(goodfit_ind,4,4);
cellEl = lbub_fits(goodfit_ind,5,4);
fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
nCells = length(goodfit_ind);
fprintf(['# goodfit cells = ' num2str(nCells) '\n'])

% input stimulus location based on experimental choice
stimEl = double(input.gratingElevationDeg);
stimAz = double(input.gratingAzimuthDeg);
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])

% calculate cell distances
fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz-stimAz).^2+(cellEl-stimEl).^2);

%% now select window to extract response
fprintf('\nExamine average cell dF/F timecourse to select response window\n')
respWindow = int64(4:12); % this will be the window to define dF/F response (offset by nITI + nDelay)

tt = (1-nOff:nOn)*(1000./input.frameImagingRateMs);
avTC = squeeze(mean(mean(tc_dfof(:,goodfit_ind,:),2),3)); % average all cells,trials
% plots all trials for one condition + one cell at a time
figure(2);clf;
for i = 1:length(goodfit_ind)
    plot(tt', squeeze(mean(tc_dfof(:,goodfit_ind(i),:),3)), 'LineWidth',1)
    hold on
end
plot(tt', avTC, 'LineWidth', 5)
title('Average cell timecourses')
ylim([-0.05 0.6])
vline(0)
vline(tt(nOff+respWindow(1)))
vline(tt(nOff+respWindow(end)))
h1=vfill([tt(nOff+respWindow(1)),tt(nOff+respWindow(end))],'gray','FaceAlpha',0.5);
uistack(h1,'bottom')

% examine individual conditions by cells
% figure(7);clf;
% for i = 1:length(goodfit_ind)
%     if ~(bIndex(i)==1)
%         continue
%     end
%     for j = 1:nSize
%         for k=1:nCon
%             iSiz = Ind_struct(j).all_trials;
%             iCon = find(conTrials == cons(k));
%             ind = intersect(iSiz, iCon);
%             
%             for l=1:length(ind)
%                 plot(tt',squeeze(tc_dfof(:,goodfit_ind(i),ind(l))))
%                 hold on
%             end
%             plot(tt',squeeze(mean(tc_dfof(:,goodfit_ind(i),ind),3)),'LineWidth',4)
%             title(['Cell #' num2str(goodfit_ind(i)) ', size ' num2str(szs(j)) ', con ' num2str(cons(k))])
%             ylim([-0.05 0.25])
%             vline(0)
%             vline(tt(nOff+window(1)))
%             vline(tt(nOff+window(end)))
%             h1=vfill([tt(nOff+window(1)),tt(nOff+window(end))],'gray','FaceAlpha',0.5);
%             pause
%             clf
%         end
%     end
% end

%% extract size tuning responses
override = 0; % over-ride for creating new data (1=override, 0=skip if exist)

filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
if exist(filename, 'file') && ~override
    fprintf('Found sizeTuneData.mat, loading previous data...\n')
    load(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    fprintf('Loaded sizeTune, sizeMean, sizeSEM, cellDists\n')
else
    fprintf('Creating new size-tuning response data...\n')
    
    % build size tuning data (with Mean and SEM) for each bin of goodfit cells
    sizeTune = cell(nSize,nCon,nCells);
    sizeMean = zeros(nSize,nCon,nCells);
    sizeSEM = sizeMean;
    
    % calculate size tuning data for all goodfit cells
    for i = 1:nCells
        iCell = goodfit_ind(i);
        
        for iCon = 1:nCon
            for iSize = 1:nSize
                ind_all = intersect(Ind_struct(iSize).all_trials,conInds{iCon});
                %stimOff = mean(mean(tc_dfof((nOff/2):nOff,iCell,ind_all),3),1);
                
                % take mean dF/F during stimOn at given ind_all
                stimOn = mean(tc_dfof(nOff+(respWindow),iCell,ind_all),1);
                sd = std(stimOn);
                
                % this cell matrix stores all stimOn for all trials in ind_all
                % dims (szs, cons, cells)
                sizeTune{iSize, iCon, i} = squeeze(stimOn);
                % also take mean and SEM
                sizeMean(iSize,iCon,i) = mean(stimOn,3);
                sizeSEM(iSize,iCon,i) = sd/sqrt(length(ind_all));
            end
        end
    end
    
    % save as sizeTuneData.mat
    save(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
    fprintf('Saved sizeTuneData.mat\n')
end

%% Fit size tuning curves
override = override_all; % over-ride for creating new data (1=override, 0=skip if exist)

% runs through all cells at each contrast condition
% calls Fit_SizeTuneSmooth_KM script which outputs fit structure

% First model: single sigmoid (3 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (6 fit params)
% use peak of fit as preferred size, and measure suppression index

verbose = 0 % verbose readout on fitting results and pause after each cell

szRng = linspace(0,max(szs));

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
%logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))

% store # trials at each size for this con con (same for all cells)
nTr = zeros(nSize,nCon);
for iSz = 1:nSize
    for iCon = 1:nCon
        nTr(iSz,iCon) = length(sizeTune{iSz,iCon,1});
    end
end

cd('K:\Code')
opts = optimset('Display','off');
    
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
if exist(filename, 'file') && ~override
    fprintf('Found sizeFitResults_SP.mat, loading previous results...\n')
    load(filename, 'sizeFits')
    fprintf('Loaded sizeFits struct\n')
else
    fprintf('Creating new size-tuning curve fit data...\n')
    
    % initialize sizeFits struct with last cell, last con
    nPts = floor(mean(nTr(:,nCon)));
    dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
    szs0 = zeros(1,nPts);%[0.1]; % and size
    for iSz = 1:nSize
        nPts = nPts + nTr(iSz,nCon);
        dumdum = [dumdum sizeTune{iSz,nCon,nCells}'];
        szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz,nCon))];
    end
    % max of each size mean for use in initial guesses
    maxMean = max(sizeMean(:,nCon,nCells));
    
    PLOTIT_FIT = 0;
    SAVEALLDATA = 1;
    Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, no plot
    %eval(['sizeFits(iCell,iCon)',' = s;']);
    
    f = fieldnames(s)';
    f{2,1} = {};
    sizeFits=struct(f{:});
    sizeFits(nCells,nCon) = s;
    
    fprintf('Begin fitting size-tuning curves at all cells, all contrasts...')
    for iCell = 1:nCells
        fprintf(['\nCell# ' num2str(iCell) '/' num2str(nCells) ', (RF-stim dist: ' num2str(cellDists(iCell)) ' deg)'])
        
        for iCon = 1:nCon
            fprintf('.')
            nPts = floor(mean(nTr(:,iCon)));
            dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
            szs0 = zeros(1,nPts);%[0.1]; % and size
            for iSz = 1:nSize
                nPts = nPts + nTr(iSz,iCon);
                dumdum = [dumdum sizeTune{iSz,iCon,iCell}'];
                szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz,iCon))];
            end
            
            % max of each size mean for use in initial guesses
            maxMean = max(sizeMean(:,iCon,iCell));
            
            PLOTIT_FIT = 0;
            SAVEALLDATA = 1;
            Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, no plot
            %eval(['sizeFits(iCell,iCon)',' = s;']);
            sizeFits(iCell,iCon)=s;
            
            if verbose
                fprintf(['\nContrast: ' num2str(cons(iCon)) ' (SStot: ' num2str(s.SStot) ')\n'])
                fprintf([num2str(nPts) ' points\n'])
                
                fprintf('Fit 1 readout:\n')
                fprintf(['Ex amp: ' num2str(s.fit1.c1(1)) ', steepness: ' num2str(s.fit1.c1(2)) ', center: ' num2str(s.fit1.c1(3))])
                fprintf(['\nFit 1 R-sq:' num2str(s.Rsq1) '\n'])
                fprintf(['Fit 1 OF:' num2str(s.fit1.OF1) ', SSE:' num2str(s.SSE1) ', Pen:' num2str(s.fit1.OF1-s.SSE1) '\n'])
                fprintf(['Pref size: ' num2str(s.prefSize1) ' (dF/F: ' num2str(0.9*s.maxResp1) ')\n'])
                fprintf(['Suppression Index: ' num2str(s.suppInd1) ' (no suppression in single sigmoid model)\n'])
                
                fprintf('Fit 2 readout:\n')
                fprintf(['Ex amp: ' num2str(s.fit2.c2(1)) ', steepness: ' num2str(s.fit2.c2(2)+s.fit2.c2(5)) ', center: ' num2str(s.fit2.c2(3))])
                fprintf(['\nInh amp: ' num2str(s.fit2.c2(4)) ', steepness: ' num2str(s.fit2.c2(5)) ', center: ' num2str(s.fit2.c2(3)+s.fit2.c2(6))])
                fprintf(['\nFit 2 R-sq: ' num2str(s.Rsq2) '\n'])
                fprintf(['Fit 2 OF:' num2str(s.fit2.OF2) ', SSE:' num2str(s.SSE2) ', Pen:' num2str(s.fit2.OF2-s.SSE2) '\n'])
                fprintf(['Pref size: ' num2str(s.prefSize2) ' (dF/F: ' num2str(s.maxResp2) ')\n'])
                fprintf(['Suppression Index: ' num2str(s.suppInd2) '\n'])
            end
        end
        
        if verbose
            fprintf(['\nF-tests for cell ' num2str(iCell) ', RF-stim dist: ' num2str(cellDists(iCell)) ' deg:\nF-scores: ' num2str([sizeFits(iCell,:).Fscore]) '\nChoose model 2? ' num2str([sizeFits(iCell,:).Ftest]) '\n'])
            pause
        end
    end
    
    % Save fit results as sizeFitResults.mat
    % save sizeFitResults.mat, with sizeFits struct
    save(filename, 'sizeFits')
    fprintf('\nSaved sizeFitResults_SP.mat\n')
end

%% scroll through individual cells fit results
if 0
    %%
    fig=figure;
    ax=axes('Parent',fig);
    i=1;
    chosen = 1:nCells; %[8 50 55 56];
    %[31 41 45 46 52 64 67 71 72 73 75 77 79 83 89]
    chosen = [29 45];
    while i < nCells+1
        if ~sum(i==chosen) %~sum(i==find(cellDists<=10)) %~sum(i==chosen)
            i = i+1;
            continue
        end
        fprintf(['Cell# ' num2str(i) '/' num2str(nCells) ', (RF-stim dist: ' num2str(cellDists(i)) ' deg)\n'])
        
        for iCon = 1:nCon
            subplot(2,nCon+2,iCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(sizeFits(i,iCon).szs0,sizeFits(i,iCon).data,'.b')
            plot(szRng,sizeFits(i,iCon).fitout1,'-')
            plot(szRng,logfit1(sizeFits(i,iCon).fit1.x0,szRng),'g--')
            hold off
            ylim([min([-0.5*sizeFits(i,iCon).maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([sizeFits(i,iCon).maxResp2 max(sizeMean(:,iCon,i))])])
            %title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits(i,iCon).Rsq1) sizeFits(i,iCon).Fstr1])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) sizeFits(i,iCon).Fstr1])
            xlabel('Stimulus Size (deg)')
            ylabel('dF/F')
            legend('mean','data','fit','Location','best')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            subplot(2,nCon+2,iCon+nCon+2)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(sizeFits(i,iCon).szs0,sizeFits(i,iCon).data,'.b')
            plot(szRng,sizeFits(i,iCon).fitout2,'-')
            plot(szRng,logfit2(sizeFits(i,iCon).fit2.x0,szRng),'g--')
            hold off
            ylim([min([-0.5*sizeFits(i,iCon).maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([sizeFits(i,iCon).maxResp2 max(sizeMean(:,iCon,i))])])
            %title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits(i,iCon).Rsq2) sizeFits(i,iCon).Fstr2])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) sizeFits(i,iCon).Fstr2])
            xlabel('Stimulus Size (deg)')
            ylabel('dF/F')
            legend('mean','data','fit','Location','best')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        
        % dF/F peak resp vs contrast - for raw data + model1 + model2
        subplot(2,nCon+2,nCon+1)
        plot(cons,max(sizeMean(:,:,i),[],1),'o')
        hold on
        plot(cons,[sizeFits(i,:).maxResp1],'^')
        plot(cons,[sizeFits(i,:).maxResp2],'s')
        hold off
        title('Peak response vs Contrast')
        xlabel('Contrast')
        ylabel('dF/F')
        legend('raw','model1','model2','Location','best')
        
        % RF-stim map
        subplot(2,nCon+2,nCon+2);cla;
        ellipse(lbub_fits(goodfit_ind(i),2,4), lbub_fits(goodfit_ind(i),3,4), 0, lbub_fits(goodfit_ind(i),4,4), lbub_fits(goodfit_ind(i),5,4),'g');
        hold on
        ellipse(min(szs)/2, min(szs)/2, 0, stimAz, stimEl, 'b');
        ellipse(max(szs)/2, max(szs)/2, 0, stimAz, stimEl, 'r');
        text(lbub_fits(goodfit_ind(i),4,4), lbub_fits(goodfit_ind(i),5,4), num2str(cellDists(i),3), 'HorizontalAlignment', 'center')
        hold off
        xlim([stimAz-50 stimAz+50])
        ylim([stimEl-50 stimEl+50])
        title('RF and stimulus map')
        xlabel('Az (deg)')
        ylabel('El (deg)')
        legend(['Cell ' num2str(i) ' - 1 sigma'], 'Min stim', 'Max stim','Location','best')
        
        % prefSize vs con - both models
        subplot(2,nCon+2,2*nCon+3)
        plot(cons,[sizeFits(i,:).prefSize1],'^')
        hold on
        plot(cons,[sizeFits(i,:).prefSize2],'s')
        plot(cons,[sizeFits(i,:).prefSize],'k-')
        hold off
        title('Pref. Size vs Contrast')
        xlabel('Contrast')
        ylabel('Pref Size (deg)')
        ylim([0 max(szs)+1])
        legend('model1','model2','Location','best')
        
        % SI vs con - both models
        subplot(2,nCon+2,2*nCon+4)
        plot(cons,[sizeFits(i,:).suppInd1],'^')
        hold on
        plot(cons,[sizeFits(i,:).suppInd2],'s')
        plot(cons,[sizeFits(i,:).suppInd],'k-')
        hold off
        title('Supp. Index vs Contrast')
        xlabel('Contrast')
        ylabel('SI')
        legend('model1','model2','Location','best')
        
        set(gcf, 'Position', [0 300 1600 500]);
        
        was_a_key = waitforbuttonpress;
        if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
            i = i - 1;
            if ~i; i=1; end;
        elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'x')
            break
        else
            i = i + 1;
        end
    end
    close(fig)
end

%% bootstrap
% bootstrap by shuffling data of highest-contrast conditions
close all;
override = override_all;

filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Fit_struct.mat']);
if exist(filename, 'file') && ~override
    fprintf('Found Fit_struct.mat, loading previous results...\n')
    load(filename, 'Fit_struct')
    fprintf('Loaded Fit_struct\n')
    Nshuf = length(Fit_struct(1).Shuf);
    fprintf(['Nshuf = ' num2str(Nshuf) '\n'])
else
    fprintf('Creating new bootstrapping data...\n')
    
    Fit_struct = [];
    
    Nshuf = 500;
    fprintf(['Nshuf = ' num2str(Nshuf) '\n'])
    
    % store # trials at each size and highest con (same for all cells)
    fprintf(['Sizes: ' num2str(szs) '\n# trials: ' num2str(nTr(:,nCon)')])
    shuf_ind = cell(nSize,1);
    
    cd('K:\Code')
    opts = optimset('Display','off');
    
    fprintf('\nBegin shuffling...\n')
    figure;
    for count_shuf = 0:Nshuf
        fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
        for iSz = 1:nSize
            if count_shuf > 0
                shuf_ind{iSz} = randsample(nTr(iSz,nCon),nTr(iSz,nCon),1); % resample with replacement
            else
                shuf_ind{iSz} = 1:nTr(iSz,nCon); % shuf_count==0 use all trials
            end
        end
        
        ifig = 1;
        start = 1;
        for iCell = 1:nCells
            %fprintf(num2str(iCell));
            % recreate sizeTune data with resampled indices
            nPts = floor(mean(nTr(:,nCon)));
            dumdum = zeros(1,nPts);%[0]; % define zero point for dF/F
            szs0 = zeros(1,nPts);%[0.1]; % and size
            for iSz = 1:nSize
                nPts = nPts + nTr(iSz,nCon);
                dum = sizeTune{iSz,nCon,iCell}';
                dumdum = [dumdum dum(shuf_ind{iSz})];
                szs0 = [szs0 szs(iSz)*ones(1,nTr(iSz,nCon))];
            end
            
            % max of each size mean for use in initial guesses
            maxMean = max(sizeMean(:,nCon,iCell));
            
            if count_shuf == 0
                PLOTIT_FIT = 1;
                SAVEALLDATA = 1;
                Fit_SizeTuneSmooth_KM % call fit script, returns fit structure s, saves plots to kevin analysis folder
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_SizeTuneSmooth_KM
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end
        if count_shuf == 0
            set(gcf, 'Position', [0 0 800 1000]);
            fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeTuneFits' num2str(ifig) '.pdf']);
            print(fn_out,'-dpdf')
        end
    end
    fprintf('\nShuffling done, saving fit results\n')
    
    save(filename, 'Fit_struct')
end

%% assess fits
% extract values for prefSize (use for decision), Ftest
% also prefSize(1/2), suppInd(1,2), fit1.c1/OF1, fit2.c2/OF2, Rsq12

fprintf('Assessing goodness of fit\n')
if Nshuf>1
    fprintf('Reading in variables of interest\n')
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)
            eval('tmp = Fit_struct(iCell).True.s_.prefSize;');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.prefSize1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.prefSize2];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.suppInd2];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Fscore];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Ftest];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.maxResp1];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.maxResp2];');
            
            % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
            fit_true_vec(iCell,:) = tmp;
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval('tmp = Fit_struct(iCell).Shuf(count_shuf).s_.prefSize;');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.prefSize1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.prefSize2];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.suppInd2];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Fscore];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Ftest];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.maxResp1];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.maxResp2];');
                
                % prefSize PS1 PS2 suppInd SI1 SI2 Fscore Ftest
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end
    
    %% plot cells with size tuning curve and shuffle results
    %chosen=[44 54]; %[31 41 45 46 52 64 67 71 72 73 75 77 79 83 89];
    %chosen = goodfit_ind_size;
    chosen = [];
    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
    ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
    for iCell = 1:nCells
        if sum(iCell==chosen)
            s = Fit_struct(iCell).True.s_;
            figure(1);clf;
            subplot(3,3,1)
            errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
            hold on
            plot(s.szs0,s.data,'.b')
            plot(szRng,s.fitout1,'-r')
            plot(szRng,s.fitout2,'-g')
            hold off
            ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
            title(['Cell #' num2str(iCell) ' Size Tuning @Con 0.8 (Ftest=' num2str(fit_true_vec(iCell,8)) ')']);
            xlabel('Stimulus size (deg)')
            ylabel('dF/F')
        end
        
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3) = mean(i); %mean
            lbub_fits(iCell,count2,5) = std(i); %stdev
            
            if sum(iCell==chosen)
                switch count2
                    case 1
                        subplot(3,3,4)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                        ylabel('count')
                    case 2
                        subplot(3,3,5)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize1 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 3
                        subplot(3,3,6)
                        histogram(i,0:max(szs)+1)
                        xlim([0 100]);
                        line([0.5*mean(i) 0.5*mean(i)], ylim, 'color','red')
                        line([2*mean(i) 2*mean(i)], ylim, 'color','red')
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['PrefSize2 shuffles cell #' num2str(iCell)])
                        xlabel('Size (deg)')
                    case 4
                        subplot(3,3,7)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd shuffles cell #' num2str(iCell)])
                        xlabel('SI')
                        ylabel('count')
                    case 5
                        subplot(3,3,8)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd1 shuffles cell #' num2str(iCell)])
                        xlabel('SI1')
                    case 6
                        subplot(3,3,9)
                        histogram(i,0:0.01:2)
                        xlim([0 1])
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        title(['SuppInd2 shuffles cell #' num2str(iCell)])
                        xlabel('SI2')
                    case 7
                        Fcrit = finv(0.95,nPts-6,nPts-3);
                        subplot(3,3,2)
                        histogram(i,50)
                        xlim([0 max(i)]);
                        line([mean(i) mean(i)], ylim)
                        line([i(ind_shuf_lb) i(ind_shuf_lb)], ylim)
                        line([i(ind_shuf_ub) i(ind_shuf_ub)], ylim)
                        line([Fcrit Fcrit], ylim, 'LineWidth', 3, 'color', 'red')
                        title(['Fscore shuffles cell #' num2str(iCell)])
                        xlabel('Fscore')
                        ylabel('count')
                    case 8
                        subplot(3,3,3)
                        histogram(i,[-0.5 0.5 1.5])
                        xlim([-0.5 1.5])
                        ylim([0 Nshuf])
                        title(['Ftest shuffles cell #' num2str(iCell)])
                        xlabel('Ftest')
                        pause
                end
            end
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:); % true (no shuffle)
    end
end

%% determine good fits
% first sort cells based on Ftest same as True fit in >50% of shuffles
% then use that model's prefSize confidence interval (e.g. prefSize1 or 2)
% check bounds of confidence interval are within 1 octave of "True" fit
% e.g. lower(2.5th)>0.5*prefSizeTrue and upper(97.5th)<2*prefSizeTrue
goodfit_ind_size = [];
for iCell = 1:nCells
    Ftest = lbub_fits(iCell,8,4); %8=Ftest, 4=True fit
    Ftestshuf = lbub_fits(iCell,8,3); %8=Ftest, 3=mean of shuffles
    lOct = 0.5*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, lower octave
    hOct = 2*lbub_fits(iCell,1,4); %1=prefSize, 4=True fit, upper octave
    switch Ftest
        case 0 % model1
            if Ftestshuf<0.5 %model1 >50% of shuffles
                if (lbub_fits(iCell,2,1)>lOct) && (lbub_fits(iCell,2,2)<hOct) %2=prefSize1, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
        case 1 % model2
            if Ftestshuf>0.5 % model2 >50% of shuffles
                if (lbub_fits(iCell,3,1)>lOct) && (lbub_fits(iCell,3,2)<hOct) %3=prefSize2, 1=lb/2=ub
                    goodfit_ind_size = [goodfit_ind_size iCell];
                end
            end
    end
end

% is model1 + is model2
ism1 = find(~lbub_fits(goodfit_ind_size,8,4));
ism2 = find(lbub_fits(goodfit_ind_size,8,4));

fprintf(['#Good cells = ' num2str(length(goodfit_ind_size)) '\nModel 1: ' num2str(length(ism1)) ', Model 2: ' num2str(length(ism2)) '\nSaving good fits\n'])

fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
save(fn_out, 'lbub_fits', 'goodfit_ind_size')

%% BOOTSTRAP plots
% ideas:
% prefSize + SI histograms (separate model1/2?)
% prefSize vs RF-distance
% SI vs RF-distance
% maxResp vs RF-distance
% maxResp vs prefSize? 
% ?

% is model1 + is model2
ism1 = find(~lbub_fits(goodfit_ind_size,8,4));
ism2 = find(lbub_fits(goodfit_ind_size,8,4));

cellDists_m1 = cellDists(goodfit_ind_size(ism1));
cellDists_m2 = cellDists(goodfit_ind_size(ism2));
prefSize_m1 = lbub_fits(goodfit_ind_size(ism1),1,4);
prefSize_m2 = lbub_fits(goodfit_ind_size(ism2),1,4);
suppInd_m1 = lbub_fits(goodfit_ind_size(ism1),4,4);
suppInd_m2 = lbub_fits(goodfit_ind_size(ism2),4,4);

% examine Ftest proportions across all cells
Ftest_allcells = lbub_fits(:,8,3); % 8=Ftest, 3=mean
figure(2);clf;
subplot(2,1,1)
histogram(Ftest_allcells,0:0.02:1)
title('Ftest proportion for all cells')
xlabel('Fraction of shuffles passing Ftest')
ylabel('Count')
subplot(2,1,2)
h1 = histogram(lbub_fits(goodfit_ind_size(ism1),8,3),0:0.02:1);
hold on
h2 = histogram(lbub_fits(goodfit_ind_size(ism2),8,3),0:0.02:1);
hold off
title('Ftest proportion for goodfit cells')
xlabel('Fraction of shuffles passing Ftest')
ylabel('Count')
legend('Model1','Model2','location','best')

% prefSize + SI
figure(3);clf;
subplot(2,1,1)
histogram(prefSize_m1,0:5:90)
hold on
histogram(prefSize_m2,0:5:90)
hold off
title('Preferred Size histogram')
xlabel('Pref Size (deg)')
ylabel('Count')
legend('Model1','Model2','location','best')
subplot(2,1,2)
histogram(suppInd_m1,0:0.05:2)
hold on
histogram(suppInd_m2,0:0.05:2)
hold off
title('Suppression Index histogram')
xlabel('SI')
ylabel('Count')
legend('Model1','Model2','location','best')

% prefSize vs RF-distance
figure(4);clf;
plot(cellDists_m1,prefSize_m1,'o')
hold on
plot(cellDists_m2,prefSize_m2,'o')
plot([0 max(szs)],[0 max(szs)],'--')
coefs1 = polyfit(cellDists_m1,prefSize_m1,1);
coefs2 = polyfit(cellDists_m2,prefSize_m2,1);
x1 = min(cellDists_m1):max(cellDists_m1);
x2 = min(cellDists_m2):max(cellDists_m2);
plot(x1,polyval(coefs1,x1),'b')
plot(x2,polyval(coefs2,x2),'r')
r1 = corrcoef(cellDists_m1,prefSize_m1);
r2 = corrcoef(cellDists_m2,prefSize_m2);
text(max(x1)+10,polyval(coefs1,max(x1)+10)+20,['y_1=' num2str(coefs1(1)) 'x+' num2str(coefs1(2))])
text(max(x2)+10,polyval(coefs2,max(x2)+10),['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x1)+10,polyval(coefs1,max(x1)+10)+15,['R^2_1=' num2str(r1(2))])
text(max(x2)+10,polyval(coefs2,max(x2)+10)-5,['R^2_2=' num2str(r2(2))])
hold off
title('Preferred Size vs RF-Stim distance')
xlabel('RF-stim dist (deg)')
ylabel('PrefSize (deg)')
legend('Model1','Model2','location','best')

% SI vs RF-distance
figure(5);clf;
plot(cellDists_m1,suppInd_m1,'o')
hold on
plot(cellDists_m2,suppInd_m2,'o')
coefs2 = polyfit(cellDists_m2,suppInd_m2,1);
x2 = min(cellDists_m2):max(cellDists_m2);
plot(x2,polyval(coefs2,x2),'r')
r2 = corrcoef(cellDists_m2,suppInd_m2);
text(max(x2)-5,polyval(coefs2,max(x2)-5)+0.2,['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x2)-5,polyval(coefs2,max(x2)-5)+0.1,['R^2_2=' num2str(r2(2))])
hold off
title('Suppression Index vs RF-Stim distance')
xlabel('RF-stim dist (deg)')
ylabel('SI')
legend('Model1','Model2','location','best')

% SI vs prefSize
figure(6);clf;
plot(prefSize_m1,suppInd_m1,'o')
hold on
plot(prefSize_m2,suppInd_m2,'o')
coefs2 = polyfit(prefSize_m2,suppInd_m2,1);
x2 = min(prefSize_m2):max(prefSize_m2);
plot(x2,polyval(coefs2,x2),'r')
r2 = corrcoef(prefSize_m2,suppInd_m2);
text(max(x2)-10,polyval(coefs2,max(x2)-10)+0.2,['y_2=' num2str(coefs2(1)) 'x+' num2str(coefs2(2))])
text(max(x2)-10,polyval(coefs2,max(x2)-10)+0.1,['R^2_2=' num2str(r2(2))])
hold off
title('Suppression Index vs Preferred Size')
xlabel('PrefSize (deg)')
ylabel('SI')
legend('Model1','Model2','location','best')

% to print
% set(gcf, 'Position', [0 0 800 1000]);
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RFs.pdf']);
% print(fn_out,'-dpdf')

%% ALL CONS FITTING PLOTS examine parameters
% these plots use data from sizeFits struct

% CHANGE: no need to bin anymore, just use cells <10deg
% compare distance to sort cells into distance ranges
fprintf('Sorting cells by RF center distance\n')
bIndex = 0*cellDists;
nBin = 3;
for i=1:nCells
    if ~sum(goodfit_ind_size == i)
        continue
    end
    if cellDists(i)<10
        bIndex(i) = 1;
    elseif cellDists(i)<15
        bIndex(i) = 2;
    else
        bIndex(i) = 3;
    end
end
nBinCells = [sum(bIndex==1) sum(bIndex==2) sum(bIndex==3)];
fprintf([num2str(nBinCells) ' cells selected\n'])

binStrs = ["0-10","10-15","15+"];
cBins = categorical({'0-10','10-15','15+'});
cBins = reordercats(cBins,{'0-10','10-15','15+'});
% binStrs = ["0-6","6-12","12+"];
% cBins = categorical({'0-6','6-12','12+'});
% cBins = reordercats(cBins,{'0-6','6-12','12+'});

% histogram of # cells in each bin
figure(1);clf;
bar(cBins,nBinCells)
title('Cell Counts by RF Distance Bin')
xlabel('Distance Bin')
ylabel('# cells')

% find cells in each bin that are also good-fits
bin1inds = find(bIndex == 1);
bin2inds = find(bIndex == 2);
bin3inds = find(bIndex == 3);


% plot all size tuning curves from each bin
for k = 1:nBin
    binCells = find(bIndex==k);
    [n, n2] = subplotn(nBinCells(k));
    figure(k+1);clf;
    suptitle(['Size Tuning Curves for cells within ' char(binStrs(k)) ' deg of stim'])
    for i = 1:nBinCells(k)
        iCell = binCells(i);
        subplot(n,n2,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean(:,iCon,iCell),sizeSEM(:,iCon,iCell))
            hold on
        end
        ylim([-1.2*(max(max(sizeSEM(:,:,iCell)))) 1.2*(max(max(sizeMean(:,:,iCell)))+max(max(sizeSEM(:,:,iCell))))])
        title(['Cell #: ' num2str(iCell) ', RFdist: ' num2str(cellDists(iCell))])
        hold off
    end
    % Construct a Legend with the data from the sub-plots
    hL = legend(num2str(cons'));
    % Programatically move the Legend
    newPosition = [0.5 0.07 0.2 0.2];
    newUnits = 'normalized';
    set(hL, 'Position', newPosition, 'Units', newUnits);
end

% last make a figure of mean size tuning curve per bin
figure(5);clf;
for k = 1:nBin
    subplot(1,nBin,k)
    for iCon=1:nCon
        errorbar(szs,mean(sizeMean(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEM(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['Average Size Tuning Curve, RF dist ' char(binStrs(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([0 0.35])
    %ylim([0 0.01])
    legend(num2str(cons'))
end

prefSize = reshape([sizeFits.prefSize], [nCells, nCon]);
prefBin1 = prefSize(bin1inds,:);
prefBin2 = prefSize(bin2inds,:);
prefBin3 = prefSize(bin3inds,:);

figure(6);clf;
subplot(1,3,1)
histogram(prefBin1,50)
title(['PrefSize (RF dist ' binStrs(1) ' deg)'])
subplot(1,3,2)
histogram(prefBin2,50)
title(['PrefSize (RF dist  ' binStrs(2) '  deg)'])
subplot(1,3,3)
histogram(prefBin3,50)
title(['PrefSize (RF dist ' binStrs(3) ' deg)'])

suppInd = reshape([sizeFits.suppInd], [nCells, nCon]);
suppBin1 = suppInd(bin1inds,:);
suppBin2 = suppInd(bin2inds,:);
suppBin3 = suppInd(bin3inds,:);
%suppBin3 = reshape([sizeFits(bin3inds,:).suppInd], [nBinCells(3), nCon]);
% truncate SI to [0,1]
suppBin1(suppBin1>1)=NaN;suppBin1(suppBin1<0)=NaN;
suppBin2(suppBin2>1)=NaN;suppBin2(suppBin2<0)=NaN;
suppBin3(suppBin3>1)=NaN;suppBin3(suppBin3<0)=NaN;

figure(7);clf;
subplot(1,3,1)
histogram(suppBin1,30)
title(['SuppInd (RF dist ' binStrs(1) ' deg)'])
subplot(1,3,2)
histogram(suppBin2,30)
title(['SuppInd (RF dist ' binStrs(2) ' deg)'])
subplot(1,3,3)
histogram(suppBin3,30)
title(['SuppInd (RF dist ' binStrs(3) ' deg)'])

% allRsq1 refers to single model, allRsq2 refers to double model

Rsq1 = reshape([sizeFits.Rsq1], [nCells, nCon]);
Rsq2 = reshape([sizeFits.Rsq2], [nCells, nCon]);
figure(8);clf;
subplot(1,3,1)
histogram(Rsq1(:),30)
ylim([0 20])
title('R^2 for single model')
subplot(1,3,2)
histogram(Rsq2(:),30)
ylim([0 20])
title('R^2 for double model')
subplot(1,3,3)
histogram(Rsq2(:)-Rsq1(:),30)
ylim([0 20])
title('Paired Difference (double-single)')

% insert adjusted Rsq *********

% now examine distance vs parameters, to pick cutoff
% ke: coefs1(2), coefs2(2)+coefs2(5)
% xe: coefs1(3),coefs2(3)
% ki: coefs2(5)
% xi: coefs2(3)+coefs2(6)
Ae1 = zeros([nCells, nCon]);
ke1 = zeros([nCells, nCon]);
xe1 = zeros([nCells, nCon]);
Ae2 = zeros([nCells, nCon]);
ke2 = zeros([nCells, nCon]);
xe2 = zeros([nCells, nCon]);
Ai2 = zeros([nCells, nCon]);
ki2 = zeros([nCells, nCon]);
xi2 = zeros([nCells, nCon]);
for i=1:nCells
    for j=1:nCon
        Ae1(i,j) = sizeFits(i,j).fit1.c1(1);
        ke1(i,j) = sizeFits(i,j).fit1.c1(2);
        xe1(i,j) = sizeFits(i,j).fit1.c1(3);
        Ae2(i,j) = sizeFits(i,j).fit2.c2(1);
        ke2(i,j) = sizeFits(i,j).fit2.c2(2)+sizeFits(i,j).fit2.c2(5);
        xe2(i,j) = sizeFits(i,j).fit2.c2(3);
        Ai2(i,j) = sizeFits(i,j).fit2.c2(4);
        ki2(i,j) = sizeFits(i,j).fit2.c2(5);
        xi2(i,j) = sizeFits(i,j).fit2.c2(3)+sizeFits(i,j).fit2.c2(6);
    end
end

x = repmat(cellDists,1,nCon)
figure(9);clf;
subplot(2,3,1)
histogram(ke1,20)
% plot(x,sizeFits.coefs1(:,:,2),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - single model')
subplot(2,3,2)
histogram(ke2,20)
% plot(x,sizeFits.coefs2(:,:,2)+sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - double model')
subplot(2,3,3)
histogram(ki2,20)
% plot(x,sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_i')
title('k_i - double model')
subplot(2,3,4)
histogram(xe1,20)
% plot(x,sizeFits.coefs1(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - single model')
subplot(2,3,5)
histogram(xe2,20)
% plot(x,sizeFits.coefs2(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - double model')
subplot(2,3,6)
histogram(xi2,20)
% plot(x,sizeFits.coefs2(:,:,3)+sizeFits.coefs2(:,:,6),'o')
% xlabel('cell dist (deg)')
% ylabel('x_i')
title('x_i - double model')

%% compare size x con fit parameters between m1 and m2, +RF-dist cutoff
% examine only in bin1 (<10 deg), compare m1 and m2 conditions
Ftest = reshape([sizeFits.Ftest], [nCells, nCon]);
ism1_all = find(~Ftest);
ism2_all = find(Ftest);
RFdistcutoff = find(repmat(cellDists, 1 ,nCon)<10);
ism1_cut = intersect(ism1_all,RFdistcutoff);
ism2_cut = intersect(ism2_all,RFdistcutoff);

% prefSize + SI
figure(10);clf;
subplot(2,1,1)
histogram(prefSize(ism1_cut),0:5:90)
hold on
histogram(prefSize(ism2_cut),0:5:90)
hold off
title('Preferred Size histogram')
xlabel('Pref Size (deg)')
ylabel('Count')
legend('Model1','Model2','location','best')
subplot(2,1,2)
histogram(suppInd(ism1_cut),0:0.05:2)
hold on
histogram(suppInd(ism2_cut),0:0.05:2)
hold off
title('Suppression Index histogram')
xlabel('SI')
ylabel('Count')
legend('Model1','Model2','location','best')

figure(11);clf;
subplot(3,4,1)
histogram(Ae1(ism1_cut),20)
hold on
histogram(Ae1(ism2_cut),20)
hold off
title('A_e - single model')
subplot(3,4,2)
histogram(Ae2(ism1_cut),20)
hold on
histogram(Ae2(ism2_cut),20)
hold off
title('A_e - double model')
subplot(3,4,3)
histogram(Ai2(ism1_cut),20)
hold on
histogram(Ai2(ism2_cut),20)
hold off
title('A_i - double model')
subplot(3,4,4)
histogram(Ae2(ism1_cut)-Ai2(ism1_cut),20)
hold on
histogram(Ae2(ism2_cut)-Ai2(ism2_cut),20)
hold off
title('A_e-A_i - double model')

subplot(3,4,5)
histogram(ke1(ism1_cut),20)
hold on
histogram(ke1(ism2_cut),20)
hold off
title('k_e - single model')
subplot(3,4,6)
histogram(ke2(ism1_cut),20)
hold on
histogram(ke2(ism2_cut),20)
hold off
title('k_e - double model')
subplot(3,4,7)
histogram(ki2(ism1_cut),20)
hold on
histogram(ki2(ism2_cut),20)
hold off
title('k_i - double model')
subplot(3,4,8)
histogram(ke2(ism1_cut)-ki2(ism1_cut),20)
hold on
histogram(ke2(ism2_cut)-ki2(ism2_cut),20)
hold off
title('k_e-k_i - double model')

subplot(3,4,9)
histogram(xe1(ism1_cut),20)
hold on
histogram(xe1(ism2_cut),20)
hold off
title('x_e - single model')
subplot(3,4,10)
histogram(xe2(ism1_cut),20)
hold on
histogram(xe2(ism2_cut),20)
hold off
title('x_e - double model')
subplot(3,4,11)
histogram(xi2(ism1_cut),20)
hold on
histogram(xi2(ism2_cut),20)
hold off
title('x_i - double model')
subplot(3,4,12)
histogram(xe2(ism1_cut)-xi2(ism1_cut),20)
hold on
histogram(xe2(ism2_cut)-xi2(ism2_cut),20)
hold off
title('x_e-x_i - double model')
legend('Model1','Model2','location','best')

% now use these matrices to examine problem cells
% ultimately define criteria for model2->model1 using low Ai2 or ki2, or another parameter
% so first, examine each cell in goodfit_ind_size, note any issues

%% examine overall data after F-test selection
cCons = categorical(cons);

% plot the number of unsuppressed vs suppressed cells
figure(12);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    subplot(1,3,k)
    
    binFtest = zeros(nCon,2); % one column each for unsupp, supp cells
    if nCutCells
        for iCon = 1:nCon
            binFtest(iCon,1) = nCutCells - sum(sizeFits.Ftest(binCells,iCon));
            binFtest(iCon,2) = sum(sizeFits.Ftest(binCells,iCon));
        end
    end
    bar(cCons,binFtest)
    title(['Unsupp vs Supp cells, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('# cells')
    legend('Unsupp','Supp')
end

figure(13);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = length(binCells);
    
    % first show average PrefSize at each contrast
    subplot(2,3,k)
    prefAv = mean(sizeFits.prefSize(binCells,:),1);
    prefSEM = std(sizeFits.prefSize(binCells,:),0,1)/sqrt(nCutCells);
    bar(cCons,prefAv)
    hold on
    errorbar(cCons,prefAv,prefSEM,'.')
    hold off
    title(['Preferred Size, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Pref Size')
    ylim([0 80])
    
    % second show average SI at each contrast
    subplot(2,3,k+3)
    SIav = nanmean(sizeFits.suppInd(binCells,:),1);
    SISEM = nanstd(sizeFits.suppInd(binCells,:),0,1)/sqrt(nCutCells);
    bar(cCons,SIav)
    hold on
    errorbar(cCons,SIav,SISEM,'.')
    hold off
    title(['Suppression Index, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Supp Ind')
    ylim([0 1])
end

figure(14);clf;
for k=1:3
    binCells = find(bIndex==k);
    
    % use only suppressed cells
    for i=1:nCon
        suppCells = intersect(binCells, find(sizeFits.Ftest(:,i)));
        nSuppCells(i) = length(suppCells);
        nSuppCellsNaN(i) = nSuppCells(i) - sum(isnan(sizeFits.suppInd(suppCells,i)));
        prefAv(i) = mean(sizeFits.prefSize(suppCells,i),1);
        prefSEM(i) = std(sizeFits.prefSize(suppCells,i))/sqrt(nSuppCells(i));
        SIav(i) = nanmean(sizeFits.suppInd(suppCells,i),1);
        SISEM(i) = nanstd(sizeFits.suppInd(suppCells,i))/sqrt(nSuppCellsNaN(i));
    end
    
    % first show average PrefSize at each contrast
    subplot(2,3,k)
    bar(cCons,prefAv)
    hold on
    errorbar(cCons,prefAv,prefSEM,'.')
    hold off
    title(['Preferred Size, RF dist: ' char(binStrs(k)) ' (n=' num2str(max(nSuppCells)) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Pref Size')
    ylim([0 80])
    
    % second show average SI at each contrast
    subplot(2,3,k+3)
    bar(cCons,SIav)
    hold on
    errorbar(cCons,SIav,SISEM,'.')
    hold off
    title(['Suppression Index, RF dist: ' char(binStrs(k)) ' (n=' num2str(max(nSuppCellsNaN)) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Supp Ind')
    ylim([0 1])
end
return

%% Examine how prefSize changes with contrast
% will need to store fit params for examining a single size

% first calculate slopes, fitting with a 1st degree polynomial (line)
lineFits = zeros(nCells,2);
for i=1:nCells
    lineFits(i,:) = polyfit(cons,(sizeFits.prefSize(i,:)/(sizeFits.prefSize(i,end))),1); %allPrefSize(i,:)
end

% the plot each cell's prefSize vs con curve, within each bin
for k=1:nBin
    figure(14+k);clf;
    [n, n2] = subplotn(nBinCells(k));
    for i = 1:nBinCells(k)
        if k==1
            iCell = bin1inds(i);
        elseif k==2
            iCell = bin2inds(i);
        else
            iCell = bin3inds(i);
        end
        subplot(n,n2,i)
        plot(cons,(sizeFits.prefSize(iCell,:)/(sizeFits.prefSize(iCell,end))),'o') %allPrefSize(iCell,:)
        hold on
        plot(cons,polyval(lineFits(iCell,:),cons),'r-')
        hold off
        title(['Cell #' num2str(iCell) ', m=' num2str(lineFits(iCell,1))])
        xlabel('Contrast')
        ylabel('Pref Sz')
        ylim([0 5])
    end
end

% readout means and stdevs for each bin
for k=1:nBin
    fprintf(['Bin #' num2str(k) ', Mean slopes: ' num2str(mean(lineFits(bIndex==k,1),1)) ', stdev: ' num2str(std(lineFits(bIndex==k,1),[],1)) '\n'])
end
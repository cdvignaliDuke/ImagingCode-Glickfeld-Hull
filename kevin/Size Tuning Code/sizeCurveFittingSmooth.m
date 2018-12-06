%% fitting with smoothness penalty
% this script uses size tuning data and fits a curve, version 3
% fit size tuning curve with either one sigmoid or sum of + and - sigmoids
% comparing by sequential F-test on each condition (cell x Con)

% this uses an expanded sigmoid model with parameters for sigmoid center
% also now using a struct to store fit results

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

% start over and load data as designated below
clear all;clc;
mouse = 'i842';
date = '180330';
ImgFolder = char('003');
time = char('1647');

RetImgFolder = char('002');

nrun = size(ImgFolder,1);
run_str = catRunName(ImgFolder, nrun);

fprintf(['Size tuning curve analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

% load behavior/experimental data
fprintf(['Loading experimental conditions from size tuning runs: ' run_str '\n'])
fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(1,:) '.mat'];
load(fName);

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
%szs = unique(Stims(:,1));
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
cellAz = lbub_fits(:,4,4);
cellEl = lbub_fits(:,5,4);
fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
nCells = length(goodfit_ind);
fprintf(['# goodfit cells = ' num2str(nCells) '\n'])

% input stimulus location based on experimental choice
stimEl = double(input.gratingElevationDeg);
stimAz = double(input.gratingAzimuthDeg);
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])

% calculate cell distances
fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz(goodfit_ind)-stimAz).^2+(cellEl(goodfit_ind)-stimEl).^2);

%% now select window to extract response
fprintf('\nExamine average cell dF/F timecourse to select response window\n')
window = int64(4:13); % this will be the window to define dF/F response (offset by nITI + nDelay)

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
vline(tt(nOff+window(1)))
vline(tt(nOff+window(end)))
h1=vfill([tt(nOff+window(1)),tt(nOff+window(end))],'gray','FaceAlpha',0.5);
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
            stimOn = mean(tc_dfof(nOff+(window),iCell,ind_all),1);
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

% plot all size tuning curves from each bin
for k = 1:nBin
    binCells = find(bIndex==k);
    [n, n2] = subplotn(nBinCells(k));
    figure(k+2);clf;
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

% save as sizeTuneData.mat
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
save(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
fprintf('Saved sizeTuneData.mat\n')

% last make a figure of mean size tuning curve per bin
figure(6);clf;
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

% %% throw out top and bottom trials, examine difference
% sizeTune_toss = cell(nSize,nCon,nCells);
% sizeMean_toss = zeros(nSize,nCon,nCells);
% sizeSEM_toss = sizeMean_toss;
% 
% % calculate size tuning data for all goodfit cells
% for i = 1:nCells
%     for iCon = 1:nCon
%         for iSize = 1:nSize
%             % load data points from sizeTune
%             stimOn = sizeTune{iSize, iCon, i};
%             
%             % find inds of min and max trials, then toss
%             minInd = find(stimOn==min(stimOn),1);
%             maxInd = find(stimOn==max(stimOn),1);
%             stimOn_toss = stimOn;
%             stimOn_toss([minInd, maxInd]) = [];
%             
%             sizeTune_toss{iSize, iCon, i} = squeeze(stimOn_toss);
%             sizeMean_toss(iSize,iCon,i) = mean(stimOn_toss);
%             sizeSEM_toss(iSize,iCon,i) = std(stimOn_toss)/sqrt(length(ind_all)-2);
%         end
%     end
% end
% 
% % plot all size tuning curves from each bin
% for k = 1:nBin
%     binCells = find(bIndex==k);
%     [n, n2] = subplotn(nBinCells(k));
%     figure(k+5);clf;
%     suptitle(['Size Tuning Curves for cells within ' char(binStrs(k)) ' deg of stim'])
%     for i = 1:nBinCells(k)
%         iCell = binCells(i);
%         subplot(n,n2,i)
%         for iCon = 1:nCon
%             errorbar(szs,sizeMean_toss(:,iCon,iCell),sizeSEM_toss(:,iCon,iCell))
%             hold on
%         end
%         ylim([-1.2*(max(max(sizeSEM_toss(:,:,iCell)))) 1.2*(max(max(sizeMean_toss(:,:,iCell)))+max(max(sizeSEM_toss(:,:,iCell))))])
%         title(['Cell #: ' num2str(iCell)])
%         hold off
%     end
%     % Construct a Legend with the data from the sub-plots
%     hL = legend(num2str(cons'));
%     % Programatically move the Legend
%     newPosition = [0.5 0.07 0.2 0.2];
%     newUnits = 'normalized';
%     set(hL, 'Position', newPosition, 'Units', newUnits);
% end
% 
% % save as sizeTuneData.mat
% % filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
% % save(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
% % fprintf('Saved sizeTuneData.mat\n')
% 
% % last make a figure of mean size tuning curve per bin
% figure(9);clf;
% for k = 1:nBin
%     subplot(2,nBin,k)
%     for iCon=1:nCon
%         errorbar(szs,mean(sizeMean(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEM(:,iCon,find(bIndex == k)),3))
%         hold on
%     end
%     title(['Average curve, RF dist ' char(binStrs(k))])
%     xlabel('Size (deg)')
%     ylabel('dF/F')
%     ylim([0 0.15])
%     ylim([0 0.01])
%     legend(num2str(cons'))
%     
%     subplot(2,nBin,k+nBin)
%     for iCon=1:nCon
%         errorbar(szs,mean(sizeMean_toss(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEM_toss(:,iCon,find(bIndex == k)),3))
%         hold on
%     end
%     title(['Average curve with toss, RF dist ' char(binStrs(k))])
%     xlabel('Size (deg)')
%     ylabel('dF/F')
%     ylim([0 0.15])
%     ylim([0 0.01])
%     legend(num2str(cons'))
% end

% averages of normalized curves in each bin (divide by peak)
% show cells re-ordered by RFdistance (maybe normalize + overlay?)

%% Fit size tuning curves
%chosen = [1];
%chosen = 1:nCells;
%chosen = [];
chosen = goodfit_ind_size;

% First model: single sigmoid (3 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (6 fit params)
% use peak of fit as preferred size, and measure suppression index

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
%logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))

szRng = linspace(0.1,max(szs));

% create a struct to store variables
s = zeros(nCells, nCon);
s3 = zeros(nCells, nCon,3);
s6 = zeros(nCells, nCon,6);
sizeFits = struct('Rsq1',s,'coefs1',s3,'Rsq2',s,'coefs2',s6,'prefSize',s,'suppInd',s,'Fscore',s,'Ftest',s);

cd('K:\Code')
opts = optimset('Display','off');

verbose = 1
figure(7);clf;
for i = 1:nCells
%     if ~bIndex(i)
%         continue
%     end
    % check if passed bootstrapping
    % here
    fprintf(['\nCell# ' num2str(i) '/' num2str(nCells) ', (RF-stim dist: ' num2str(cellDists(i)) ' deg)'])
    
    for iCon = 1:nCon
        dumdum = zeros(1,8);%[0]; % define zero point for dF/F
        szs0 = zeros(1,8);%[0.1]; % and size
        nPts = 8;
        for iSz = 1:nSize
            nPts = nPts + length(sizeTune{iSz,iCon,i});
            dumdum = [dumdum sizeTune{iSz,iCon,i}'];
            szs0 = [szs0 szs(iSz)*ones(1,length(sizeTune{iSz,iCon,i}))];
        end
        SStot = sum((dumdum-mean(dumdum)).^2);
        
        if verbose
            fprintf(['\nContrast: ' num2str(cons(iCon))])
            fprintf([' (SStot: ' num2str(SStot) ')'])
            fprintf('\nEvaluating fit 1...\n')
        end
        
        % max of each size mean for use in initial guesses
        maxMean = max(sizeMean(:,iCon,i));
        
        % single sigmoid model
        guess1 = [0.5*maxMean 0.35 10]; % guesses for 3 params (Ae, ke, xe)
        lb = [0 0 0];
        ub = [inf 1 80]; %[inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        %[fitx,resnorm1] = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        %guess = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        [fitx,OF1] = fminsearchbnd(@(c)sigevalpen1(c,szs0,dumdum),guess1,lb,ub,opts);
        resnorm1 = sigeval1(fitx,szs0,dumdum);
        r2 = 1 - resnorm1/SStot; % resnorm is SSE
        sizeFits.Rsq1(i, iCon) = r2;
        sizeFits.coefs1(i,iCon,:) = fitx;
        
        fitout1 = logfit1(fitx,szRng);
        maxResp1 = max(fitout1);
        prefSize1 = szRng(find(fitout1>(0.9*maxResp1),1));
        if ~length(prefSize1)
            prefSize1=max(szs);
        end
        suppInd1 = 0;
        
        if verbose
            fprintf('Fit 1 readout:\n')
            fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)) ', center: ' num2str(fitx(3))])
            fprintf(['\nFit 1 R-sq:' num2str(r2) '\n'])
            fprintf(['Fit 1 OF:' num2str(OF1) ', SSE:' num2str(resnorm1) ', Pen:' num2str(OF1-resnorm1) '\n'])
            fprintf(['Pref size: ' num2str(prefSize1) ' (dF/F: ' num2str(0.9*maxResp1) ')\n'])
            fprintf(['Suppression Index: ' num2str(suppInd1) ' (no suppression in single sigmoid model)\n'])
            fprintf('\nEvaluating fit 2...\n')
        end
        
        % double sigmoid model
        guess2 = [maxMean 0.4 10 maxMean-mean(dumdum) 0.3 20]; % guesses for (Ae, k1, x1, Ai, k2, x2)
        lb = [0 0 0 0 0 0];
        ub = [inf 1 80 inf 1 80]; %[inf inf inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        %[fitx,resnorm2] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        %guess = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        [fitx,OF2] = fminsearchbnd(@(c)sigevalpen2(c,szs0,dumdum),guess2,lb,ub,opts);
        resnorm2 = sigeval2(fitx,szs0,dumdum);
        r2 = 1 - resnorm2/SStot; % resnorm is SSE
        sizeFits.Rsq2(i,iCon) = r2;
        sizeFits.coefs2(i,iCon,:) = fitx;
        
        fitout2 = logfit2(fitx,szRng);
        maxResp2 = max(fitout2);
        prefSize2 = szRng(find(fitout2==maxResp2,1));
        suppInd2 = 1 - fitout2(end)/maxResp2;
        
        if verbose
            fprintf('Fit 2 readout:\n')
            fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(5)) ', center: ' num2str(fitx(3))])
            fprintf(['\nInh amp: ' num2str(fitx(4)) ', steepness: ' num2str(fitx(5)) ', center: ' num2str(fitx(3)+fitx(6))])
            fprintf(['\nFit 2 R-sq: ' num2str(r2) '\n'])
            fprintf(['Fit 2 OF:' num2str(OF2) ', SSE:' num2str(resnorm2) ', Pen:' num2str(OF2-resnorm2) '\n'])
            fprintf(['Pref size: ' num2str(prefSize2) ' (dF/F: ' num2str(maxResp2) ')\n'])
            fprintf(['Suppression Index: ' num2str(suppInd2) '\n'])
        end
        
        % plot both fits, if chosen cell
        if sum(i==chosen)
            subplot(2,nCon,iCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout1,'-')
            plot(szRng,logfit1(guess1,szRng),'g--')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq1(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('mean','data','fit')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            subplot(2,nCon,iCon+nCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout2,'-')
            plot(szRng,logfit2(guess2,szRng),'g--')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq2(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('mean','data','fit')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        
        % now F-test
        if verbose
            fprintf([num2str(nPts) ' points\n'])
        end
        dfS = nPts - 3; % single sigmoid model df = #sizes - 3 model parameters
        dfD = nPts - 6; % double sigmoid model df = #sizes - 6 model parameters
        Fcrit = finv(0.95,dfD,dfS); % measure critical F at a=0.05
        sizeFits.Fscore(i,iCon) = (resnorm2 - resnorm1)/(dfD-dfS) ./ (resnorm2/dfD);
        Ftest = sizeFits.Fscore(i,iCon)>Fcrit;
        sizeFits.Ftest(i,iCon) = Ftest;
        
        % store values from selected test (0=single,1=double)
        if Ftest
            sizeFits.prefSize(i,iCon) = prefSize2;
            sizeFits.suppInd(i,iCon) = suppInd2;
        else
            sizeFits.prefSize(i,iCon) = prefSize1;
            sizeFits.suppInd(i,iCon) = suppInd1;
        end
    end
    
    if verbose
        fprintf(['\nF-tests for cell i=' num2str(i) ', RF-stim dist: ' num2str(cellDists(i)) ' deg:\nF-scores: ' num2str(sizeFits.Fscore(i,:)) '\nChoose model 2? ' num2str(sizeFits.Ftest(i,:)) '\n'])
    end
    
    if sum(i==chosen)
        pause
    end
end

% Save fit results as sizeFitResults.mat
% save sizeFitResults.mat, with sizeFits struct, bIndex
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
save(filename, 'sizeFits')
fprintf('\nSaved sizeFitResults_SP.mat\n')

return

%% Load instead of fitting
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
load(filename)
fprintf('Loaded sizeTuneData.mat\n')

filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
load(filename)
fprintf('Loaded sizeFitResults_SP.mat\n')

% WARNING: bIndex is calculated based on predefined bin ranges (right now 0-6, 6-12, 12+)
% so better to use cellDists and recalculate bins if they change

%% examine parameters
% compare distance to index cells into distance ranges
fprintf('Sorting cells by RF center distance\n')
bIndex = 0*cellDists;
nBin = 3;
for i=1:nCells
    if cellDists(i)<7
        bIndex(i) = 1;
    elseif cellDists(i)<15
        bIndex(i) = 2;
    else
        bIndex(i) = 3;
    end
end
nBinCells = [sum(bIndex==1) sum(bIndex==2) sum(bIndex==3)];
fprintf([num2str(nBinCells) ' cells selected\n'])

binStrs = ["0-7","7-15","15+"];
cBins = categorical({'0-7','7-15','15+'});
cBins = reordercats(cBins,{'0-7','7-15','15+'});
% binStrs = ["0-6","6-12","12+"];
% cBins = categorical({'0-6','6-12','12+'});
% cBins = reordercats(cBins,{'0-6','6-12','12+'});

% histogram of # cells in each bin
figure(1);clf;
bar(cBins,nBinCells)
title('Cell Counts by RF Distance Bin')
xlabel('Distance Bin')
ylabel('# cells')

bin1inds = find(bIndex == 1);
bin2inds = find(bIndex == 2);
bin3inds = find(bIndex == 3);

figure(3);clf;
subplot(1,3,1)
histogram(sizeFits.prefSize(bin1inds,:),50)
title('PrefSize (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(sizeFits.prefSize(bin2inds,:),50)
title('PrefSize (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(sizeFits.prefSize(bin3inds,:),50)
title('PrefSize (RF dist 15+ deg)')

% truncate SI to [0,1]
SI=sizeFits.suppInd;SI(SI>1)=NaN;SI(SI<0)=NaN;
sizeFits.suppInd=SI;

figure(4);clf;
subplot(1,3,1)
histogram(sizeFits.suppInd(bin1inds,:),30)
title('SuppInd (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(sizeFits.suppInd(bin2inds,:),30)
title('SuppInd (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(sizeFits.suppInd(bin3inds,:),30)
title('SuppInd (RF dist 15+ deg)')

% allRsq1 refers to single model, allRsq2 refers to double model
figure(5);clf;
subplot(1,3,1)
histogram(sizeFits.Rsq1(:),30)
ylim([0 20])
title('R^2 for single model')
subplot(1,3,2)
histogram(sizeFits.Rsq2(:),30)
ylim([0 20])
title('R^2 for double model')
subplot(1,3,3)
histogram(sizeFits.Rsq2(:)-sizeFits.Rsq1(:),30)
ylim([0 20])
title('Paired Difference (double-single)')

% insert adjusted Rsq *********

% now examine distance vs parameters, to pick cutoff
% ke: coefs1(2), coefs2(2)+coefs2(5)
% xe: coefs1(3),coefs2(3)
% ki: coefs2(5)
% xi: coefs2(3)+coefs2(6)
x = repmat(cellDists,1,nCon)
figure(6);clf;
subplot(2,3,1)
histogram(sizeFits.coefs1(:,:,2),20)
% plot(x,sizeFits.coefs1(:,:,2),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - single model')
subplot(2,3,2)
histogram(sizeFits.coefs2(:,:,2)+sizeFits.coefs2(:,:,5),20)
% plot(x,sizeFits.coefs2(:,:,2)+sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - double model')
subplot(2,3,3)
histogram(sizeFits.coefs2(:,:,5),20)
% plot(x,sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_i')
title('k_i - double model')
subplot(2,3,4)
histogram(sizeFits.coefs1(:,:,3),20)
% plot(x,sizeFits.coefs1(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - single model')
subplot(2,3,5)
histogram(sizeFits.coefs2(:,:,3),20)
% plot(x,sizeFits.coefs2(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - double model')
subplot(2,3,6)
histogram(sizeFits.coefs2(:,:,3)+sizeFits.coefs2(:,:,6),20)
% plot(x,sizeFits.coefs2(:,:,3)+sizeFits.coefs2(:,:,6),'o')
% xlabel('cell dist (deg)')
% ylabel('x_i')
title('x_i - double model')

%% examine overall data after F-test selection
cCons = categorical(cons);

% plot the number of unsuppressed vs suppressed cells
figure(7);clf;
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

figure(8);clf;
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

figure(9);clf;
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
    figure(8+k);clf;
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
%% this script uses size tuning data and fits a curve
% fit size tuning curve with either one sigmoid or sum of + and - sigmoids
% comparing by sequential F-test on each condition (cell x Con)

% second version, looking at all cells by RF distance instead of just <7 deg

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

loadFlag = 1;
if loadFlag
    % start over and load data as designated below
    clear all;clc;
    mouse = 'i908';
    date = '171229';
    ImgFolder = char('003','004','005','006');
    time = char('0923','0953','1021','1053');
    
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
end

% input stimulus location based on experimental choice
stimEl = 10;
stimAz = 20;
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])

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

%% start analyzing cell distances
fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz(goodfit_ind)-stimAz).^2+(cellEl(goodfit_ind)-stimEl).^2);

% compare distance to index cells into distance ranges
fprintf('Sorting cells by RF center distance\n')
bIndex = 0*cellDists;
nBin = 3;
for i=1:nCells
    if cellDists(i)<6
        bIndex(i) = 1;
    elseif cellDists(i)<12
        bIndex(i) = 2;
    else
        bIndex(i) = 3;
    end
end
nBinCells = [sum(bIndex==1) sum(bIndex==2) sum(bIndex==3)];
fprintf([num2str(nBinCells) ' cells selected\n'])

% binStrs = ["0-7","7-15","15+"];
% cBins = categorical({'0-7','7-15','15+'});
% cBins = reordercats(cBins,{'0-7','7-15','15+'});
binStrs = ["0-6","6-12","12+"];
cBins = categorical({'0-6','6-12','12+'});
cBins = reordercats(cBins,{'0-6','6-12','12+'});

% histogram of # cells in each bin
figure(1);clf;
bar(cBins,nBinCells)
title('Cell Counts by RF Distance Bin')
xlabel('Distance Bin')
ylabel('# cells')

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
            stimOn = mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),1);
            sd = std(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),1));
            
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
    figure(k+1);clf;
    suptitle(['Size Tuning Curves for cells within ' char(binStrs(k)) ' deg of stim'])
    for i = 1:nBinCells(k)
        subplot(n,n2,i)
        for iCon = 1:nCon
            errorbar(szs,sizeMean(:,iCon,i),sizeSEM(:,iCon,i))
            hold on
        end
        ylim([-1.2*(max(max(sizeSEM(:,:,i)))) 1.2*(max(max(sizeMean(:,:,i)))+max(max(sizeSEM(:,:,i))))])
        title(['Cell #: ' num2str(i)])
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
    ylim([0 0.15])
    legend(num2str(cons'))
end

%% Fit size tuning curves
chosen = 131; %129 131

% First model: single sigmoid (2 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (4 fit params)
% use peak of fit as preferred size, and measure suppression index

% single sigmoid fits Ae, ke=k1
logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, Ai, ki=k2
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(4))*(xdata))) - coefs(3)./(1+exp(-coefs(4)*(xdata)))

szRng = linspace(0.1,max(szs));

% create variable vectors for data
allRsq1 = zeros(nCells, nCon); % cells x contrasts, collects rsq for each fit
allRsq2 = allRsq1;
allPrefSize = allRsq1; % these two are for final Pref and SI, after F-test
allSuppInd = allRsq1;
allFscore = allRsq1; % F-score for comparing two fits at each cell x con
allFtest = allFscore; % results of F-test (0 = single fit, 1 = double fit)

opts = optimset('Display','off');

figure(6);clf;
for i = 1:nCells
    if ~bIndex(i)
        continue
    end
    fprintf(['\nCell# ' num2str(i) ', bIndex ' num2str(bIndex(i))])
    
    for iCon = 1:nCon
        fprintf(['\nContrast: ' num2str(cons(iCon))])
        
        dumdum = [0]; % define zero point for dF/F
        szs0 = [0.1]; % and size
        nPts = 1;
        for iSz = 1:nSize
            nPts = nPts + length(sizeTune{iSz,iCon,i});
            dumdum = [dumdum sizeTune{iSz,iCon,i}'];
            szs0 = [szs0 szs(iSz)*ones(1,length(sizeTune{iSz,iCon,i}))];
        end
        
        SStot = sum((dumdum-mean(dumdum)).^2);
        
        % single sigmoid model
        fprintf('\nEvaluating fit 1...\n')
        guess = [max(dumdum) 0.35]; % guesses for 2 params (Ae, ke)
        lb = [0 0];
        ub = [inf 1]; %[inf inf]; % here limit steepness to 1 (above this is way too steep)
        [fitx,resnorm1] = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        
        fprintf('Fit 1 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2))])
        r2 = 1 - resnorm1/SStot; % resnorm is SSE
        fprintf(['\nFit 1 R-sq:' num2str(r2) '\n'])
        allRsq1(i, iCon) = r2;
        %allSSE1(i, iCon) = resnorm1;
        
        fitout1 = logfit1(fitx,szRng);
        maxResp1 = max(fitout1);
        prefSize1 = szRng(find(fitout1>(0.9*maxResp1),1));
        suppInd1 = 0;
        fprintf(['Pref size: ' num2str(prefSize1) ' (dF/F: ' num2str(0.9*maxResp1) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd1) ' (no suppression in single sigmoid model)\n'])
        
        % double sigmoid model
        fprintf('\nEvaluating fit 2...\n')
        guess = [2*max(dumdum) 0.1 mean(dumdum) 0.06]; % guesses for (Ae, k1, Ai, k2)
        lb = [0 0 0 0];
        ub = [inf 1 inf 1]; %[inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        [fitx,resnorm2] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        
        fprintf('Fit 2 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(4))])
        fprintf(['\nInh amp: ' num2str(fitx(3)) ', steepness: ' num2str(fitx(4))])
        r2 = 1 - resnorm2/SStot; % resnorm is SSE
        fprintf(['\nFit 2 R-sq: ' num2str(r2) '\n'])
        allRsq2(i, iCon) = r2;
        %allSSE2(i, iCon) = resnorm2;
        
        fitout2 = logfit2(fitx,szRng);
        maxResp2 = max(fitout2);
        prefSize2 = szRng(find(fitout2==maxResp2,1));
        suppInd2 = 1 - fitout2(end)/maxResp2;
        fprintf(['Pref size: ' num2str(prefSize2) ' (dF/F: ' num2str(maxResp2) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd2) '\n'])
        
        % plot both fits, if chosen cell
        if 0%i==chosen
            subplot(2,nCon,iCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout1,'-')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(allRsq1(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
            
            subplot(2,nCon,iCon+nCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout2,'-')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(allRsq2(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
        end
        
        % now F-test
        fprintf([num2str(nPts) ' points\n'])
        dfS = nPts - 2; % single sigmoid model df = #sizes - 2 model parameters
        dfD = nPts - 4; % double sigmoid model df = #sizes - 4 model parameters
        Fcrit = finv(0.95,dfD,dfS); % measure critical F at a=0.05
        allFscore(i,iCon) = (resnorm2 - resnorm1)/(dfD-dfS) ./ (resnorm2/dfD);
        Ftest = allFscore(i,iCon)>Fcrit;
        allFtest(i,iCon) = Ftest;
        
        % store values from selected test (0=single,1=double)
        if Ftest
            allPrefSize(i,iCon) = prefSize2;
            allSuppInd(i,iCon) = suppInd2;
        else
            allPrefSize(i,iCon) = prefSize1;
            allSuppInd(i,iCon) = suppInd1;
        end
    end
    
    fprintf(['\nF-tests for cell i=' num2str(i) ', distIndex ' num2str(bIndex(i)) ':\nF-scores: ' num2str(allFscore(i,:)) '\nChoose model 2? ' num2str(allFtest(i,:)) '\n'])
    
    if 0%bIndex(i)==2
        pause
    end
end

% Save fit results as sizeFitResults.mat
% save sizeFitResults.mat, with allFtest, allPrefSize, allSuppInd, bIndex
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults.mat']);
save(filename, 'allFtest', 'allPrefSize', 'allSuppInd', 'bIndex')
fprintf('Saved sizeFitResults.mat\n')

%% Load instead of fitting
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults.mat']);
load(filename,'allFtest', 'allPrefSize', 'allSuppInd') %, 'bIndex')
fprintf('Loaded sizeFitResults.mat\n')

%% examine parameters
bin1inds = find(bIndex == 1);
bin2inds = find(bIndex == 2);
bin3inds = find(bIndex == 3);

figure(3);clf;
subplot(1,3,1)
histogram(allPrefSize(bin1inds,:),50)
title('PrefSize (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(allPrefSize(bin2inds,:),50)
title('PrefSize (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(allPrefSize(bin3inds,:),50)
title('PrefSize (RF dist 15+ deg)')

figure(4);clf;
subplot(1,3,1)
histogram(allSuppInd(bin1inds,:),30)
title('SuppInd (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(allSuppInd(bin2inds,:),30)
title('SuppInd (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(allSuppInd(bin3inds,:),30)
title('SuppInd (RF dist 15+ deg)')

% allRsq1 refers to single model, allRsq2 refers to double model
figure(5);clf;
subplot(1,3,1)
histogram(allRsq1(:),30)
ylim([0 20])
title('R^2 for single model')
subplot(1,3,2)
histogram(allRsq2(:),30)
ylim([0 20])
title('R^2 for double model')
subplot(1,3,3)
histogram(allRsq2(:)-allRsq1(:),30)
ylim([0 20])
title('Paired Difference (double-single)')

% insert adjusted Rsq *********

%% examine overall data after F-test selection
cCons = categorical(cons);

% add a plot for the number of unsuppressed vs suppressed cells
figure(6);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    subplot(1,3,k)
    
    binFtest = zeros(nCon,2); % one column each for unsupp, supp cells
    if nCutCells
        for iCon = 1:nCon
            binFtest(iCon,1) = nCutCells - sum(allFtest(binCells,iCon));
            binFtest(iCon,2) = sum(allFtest(binCells,iCon));
        end
    end
    bar(cCons,binFtest)
    title(['Unsupp vs Supp cells, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('# cells')
    legend('Unsupp','Supp')
end

figure(7);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = length(binCells);
    
    % first show average PrefSize at each contrast
    subplot(2,3,k)
    prefAv = mean(allPrefSize(binCells,:),1);
    prefSEM = std(allPrefSize(binCells,:),0,1)/sqrt(nCutCells);
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
    SIav = mean(allSuppInd(binCells,:),1);
    SISEM = std(allSuppInd(binCells,:),0,1)/sqrt(nCutCells);
    bar(cCons,SIav)
    hold on
    errorbar(cCons,SIav,SISEM,'.')
    hold off
    title(['Suppression Index, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Supp Ind')
    ylim([0 0.6])
end

%% Examine how prefSize changes with contrast
% will need to store fit params for examining a single size

% first calculate slopes, fitting with a 1st degree polynomial (line)
lineFits = zeros(nCells,2);
for i=1:nCells
    lineFits(i,:) = polyfit(cons,allPrefSize(i,:),1);
end

% the plot each cell's prefSize vs con curve, within each bin
for k=1:nBin
    figure(7+k);clf;
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
        plot(cons,allPrefSize(iCell,:),'o')
        hold on
        plot(cons,polyval(lineFits(iCell,:),cons),'r-')
        hold off
        title(['Cell #' num2str(iCell) ', m=' num2str(lineFits(iCell,1))])
        xlabel('Contrast')
        ylabel('Pref Sz')
        ylim([0 90])
    end
end

% readout means and stdevs for each bin
for k=1:nBin
    fprintf(['Bin #' num2str(k) ', Mean slopes: ' num2str(mean(lineFits(bIndex==k,1),1)) ', stdev: ' num2str(std(lineFits(bIndex==k,1),[],1)) '\n'])
end
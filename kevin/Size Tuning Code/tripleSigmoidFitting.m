%% size tuning fitting with 3 sigmoids
% this script uses size tuning data and fits a curve, version 3
% fit size tuning curve with one of:
% one + sigmoid
% sum of + and - sigmoids (2)
% sum of + - and + sigmoids (3)
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

loadFlag = 1;
if loadFlag
    % start over and load data as designated below
    clear all;clc;
    mouse = 'i926';
    date = '161206';
    ImgFolder = char('004','006');
    time = char('1532','1621');
    
    RetImgFolder = char('001');
    
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
stimEl = -15;
stimAz = 5;
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

% add code to throw out top and bottom trial

%% throw out top and bottom trials, examine difference

%% Fit size tuning curves
chosen = 131;

% First model: single sigmoid (3 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (6 fit params)
% use peak of fit as preferred size, and measure suppression index
% Third model: sum of +, - and + sigmoids (9 fit params)
% use peak of fit as preferred size, and measure suppression index

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
%logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) ...
                       - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))
% triple sigmoid fits Ae1, ke1=k1+k2+k3, xe=x1, Ai, ki=k2+k3, xi=x1+x2, Ae2, ke2=k3, xe2=x1+x2+x3
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit3 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5)+coefs(8))*(xdata-coefs(3)))) ...
                       - coefs(4)./(1+exp(-(coefs(5)+coefs(8))*(xdata-(coefs(3)+coefs(6))))) ...
                       + coefs(7)./(1+exp(-coefs(8)*(xdata-(coefs(3)+coefs(6)+coefs(9)))))

szRng = linspace(0.1,max(szs));

% create a struct to store variables
s = zeros(nCells, nCon);
s3 = zeros(nCells, nCon,3);
s6 = zeros(nCells, nCon,6);
s9 = zeros(nCells, nCon,9);
sizeFits = struct('Rsq1',s,'coefs1',s3,'Rsq2',s,'coefs2',s6,'Rsq3',s,'coefs3',s9,'prefSize',s,'suppInd',s,'Fscore1',s,'Fscore2',s,'Ftest',s);

opts = optimset('Display','off'); % turn verbose fitting off

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
        guess = [max(dumdum) 0.35 10]; % guesses for 3 params (Ae, ke, xe)
        lb = [0 0 0];
        ub = [inf 1 100]; %[inf inf]; % here limit steepness to 1 (above this is way too steep)
        [fitx,resnorm1] = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        
        fprintf('Fit 1 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)) ', center: ' num2str(fitx(3))])
        r2 = 1 - resnorm1/SStot; % resnorm is SSE
        fprintf(['\nFit 1 R-sq:' num2str(r2) '\n'])
        sizeFits.Rsq1(i, iCon) = r2;
        sizeFits.coefs1(i,iCon,:) = fitx;
        %allSSE1(i, iCon) = resnorm1;
        
        fitout1 = logfit1(fitx,szRng);
        maxResp1 = max(fitout1);
        prefSize1 = szRng(find(fitout1>(0.9*maxResp1),1));
        suppInd1 = 0;
        fprintf(['Pref size: ' num2str(prefSize1) ' (dF/F: ' num2str(0.9*maxResp1) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd1) ' (no suppression in single sigmoid model)\n'])
        
        % double sigmoid model
        fprintf('\nEvaluating fit 2...\n')
        guess = [2*max(dumdum) 0.1 10 mean(dumdum) 0.06 20]; % guesses for (Ae, k1, x1, Ai, k2 x2)
        lb = [0 0 0 0 0 0];
        ub = [inf 1 100 inf 1 100]; %[inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        [fitx,resnorm2] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        
        fprintf('Fit 2 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(5)) ', center: ' num2str(fitx(3))])
        fprintf(['\nInh amp: ' num2str(fitx(4)) ', steepness: ' num2str(fitx(5)) ', center: ' num2str(fitx(3)+fitx(6))])
        r2 = 1 - resnorm2/SStot; % resnorm is SSE
        fprintf(['\nFit 2 R-sq: ' num2str(r2) '\n'])
        sizeFits.Rsq2(i,iCon) = r2;
        sizeFits.coefs2(i,iCon,:) = fitx;
        %allSSE2(i, iCon) = resnorm2;
        
        fitout2 = logfit2(fitx,szRng);
        maxResp2 = max(fitout2);
        prefSize2 = szRng(find(fitout2==maxResp2,1));
        suppInd2 = 1 - fitout2(end)/maxResp2;
        fprintf(['Pref size: ' num2str(prefSize2) ' (dF/F: ' num2str(maxResp2) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd2) '\n'])
        
        % triple sigmoid model
        fprintf('\nEvaluating fit 3...\n')
        guess = [2*max(dumdum) 0.1 10 mean(dumdum) 0.06 20 0.5*mean(dumdum) 0.01 30]; % guesses for (Ae1 k1 x1 Ai k2 x2 Ae2 k3 x3)
        lb = [0 0 0 0 0 0 0 0 0];
        ub = [inf 1 100 inf 1 100 inf 1 100]; %[inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        [fitx,resnorm3] = lsqcurvefit(logfit3,guess,szs0,dumdum,lb,ub,opts);
        
        fprintf('Fit 3 readout:\n')
        fprintf(['Ex1 amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(5)+fitx(8)) ', center: ' num2str(fitx(3))])
        fprintf(['\nInh amp: ' num2str(fitx(4)) ', steepness: ' num2str(fitx(5)+fitx(8)) ', center: ' num2str(fitx(3)+fitx(6))])
        fprintf(['\nEx2 amp: ' num2str(fitx(7)) ', steepness: ' num2str(fitx(8)) ', center: ' num2str(fitx(3)+fitx(6)+fitx(9))])
        r2 = 1 - resnorm3/SStot; % resnorm is SSE
        fprintf(['\nFit 3 R-sq: ' num2str(r2) '\n'])
        sizeFits.Rsq3(i,iCon) = r2;
        sizeFits.coefs3(i,iCon,:) = fitx;
        %allSSE2(i, iCon) = resnorm2;
        
        fitout3 = logfit3(fitx,szRng);
        maxResp3 = max(fitout3);
        prefSize3 = szRng(find(fitout3==maxResp3,1));
        suppInd3 = 1 - fitout3(end)/maxResp3;
        fprintf(['Pref size: ' num2str(prefSize3) ' (dF/F: ' num2str(maxResp3) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd3) '\n'])
        
        % plot both fits, if chosen cell
        if i==50
            subplot(3,nCon,iCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout1,'-')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq1(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
            
            subplot(3,nCon,iCon+nCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout2,'-')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq2(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
            
            subplot(3,nCon,iCon+2*nCon)
            errorbar([0 szs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout3,'-')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq3(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
        end
        
        % now F-test
        fprintf([num2str(nPts) ' points\n'])
        dfS = nPts - 3; % single sigmoid model df = #sizes - 3 model parameters
        dfD = nPts - 6; % double sigmoid model df = #sizes - 6 model parameters
        dfT = nPts - 9; % triple sigmoid model df = #sizes - 9 model parameters
        Fcrit1 = finv(0.95,dfD,dfS); % measure critical F at a=0.05
        sizeFits.Fscore1(i,iCon) = (resnorm2 - resnorm1)/(dfD-dfS) ./ (resnorm2/dfD);
        Fcrit2 = finv(0.95,dfT,dfD); % measure critical F at a=0.05
        sizeFits.Fscore2(i,iCon) = (resnorm3 - resnorm2)/(dfT-dfD) ./ (resnorm3/dfT);
        Ftest = sizeFits.Fscore1(i,iCon)>Fcrit1; % Ftest=1 if Fscore1>Fcrit1
        Ftest = Ftest + Ftest.*(sizeFits.Fscore2(i,iCon)>Fcrit2); % Ftest=2 if Fscore2>Fcrit2 and Ftest==1
        sizeFits.Ftest(i,iCon) = Ftest;
        
        if i==50
            pause
        end
        % store values from selected test (0=single,1=double,2=triple)
        if Ftest==2
            sizeFits.prefSize(i,iCon) = prefSize3;
            sizeFits.suppInd(i,iCon) = suppInd3;
        elseif Ftest==1
            sizeFits.prefSize(i,iCon) = prefSize2;
            sizeFits.suppInd(i,iCon) = suppInd2;
        else
            sizeFits.prefSize(i,iCon) = prefSize1;
            sizeFits.suppInd(i,iCon) = suppInd1;
        end
    end
    
    fprintf(['\nF-tests for cell i=' num2str(i) ', distIndex ' num2str(bIndex(i)) ':\nF-tests: ' num2str(sizeFits.Ftest(i,:)) '\n'])
    
    if i==50
        pause
    end
end

% Save fit results as sizeFitResults.mat
% save sizeFitResults.mat, with sizeFits struct, bIndex
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_v3.mat']);
save(filename, 'sizeFits', 'bIndex')
% save(filename, 'allFtest', 'allPrefSize', 'allSuppInd', 'bIndex')
fprintf('Saved sizeFitResults.mat\n')

%% Load instead of fitting
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
load(filename)
fprintf('Loaded sizeTuneData.mat\n')

filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults.mat']);
load(filename) %'allFtest', 'allPrefSize', 'allSuppInd') %, 'bIndex')
fprintf('Loaded sizeFitResults.mat\n')

% WARNING: bIndex is calculated based on predefined bin ranges (right now 0-6, 6-12, 12+)
% so better to use cellDists and recalculate bins if they change

%% examine parameters
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
    SIav = mean(sizeFits.suppInd(binCells,:),1);
    SISEM = std(sizeFits.suppInd(binCells,:),0,1)/sqrt(nCutCells);
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
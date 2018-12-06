%% Load and Combine Size Tuning Data
% script by KM, updated: 3/26/2018

%% initialize, select data
% currently only choose data with 8 size x 4 con, fitting done
% will add functionality later for 8 size x 3 con (incorporate Valerie's data)

clear all;clc;

fprintf('Size tuning multi-dataset analysis - by KM, Glickfeld Lab\nSelected data:\n')
%             %mouse  date      ret_run sz_run;
dataSelect = [%"i901", "171204", "002", "003";
              "i901" "171226", "002", "004-007";
              %"i901" "180316", "002", "003"; % PM, not great
              %"i904" "171227", "002", "003";
              %"i906" "171216", "002", "003-004"; %-005 if fix frame issue
              "i906" "171228", "002", "003-006";
              "i906" "171230", "002", "003-006";
              % "i908" "171207", "003", "004-005"; % bad data
              "i908" "171229", "002", "003-006";
              "i908" "180108", "002", "003-004";
              "i908" "180208", "002", "003-004";
              "i927" "180227", "002", "003-004";
              "i927" "180308", "002", "003-004";
              "i927" "180312", "002", "003-004";
              % "i842" "180326", "003", "004"; % GCaMP6s
              ]
          
%                   %mouse  date      ret_run sz_run    stimEl stimAz;
% dataSelect3con = [%"i901", "171204", "002", "003" 0 0; % not 0 0
%                   %"i901" "171226", "002", "004-007" 10 20;
%                   %"i904" "171227", "002", "003" 0 15;
%                   %"i906" "171216", "002", "003-004" 10 10; %-005? if fix behav data frames
%                   %"i906" "171228", "002", "003-006" 10 5;
%                   "i906" "171230", "002", "003-006" 15 5;
%                   % "i908" "171207", "002", "004-005" 15 10;
%                   "i908" "171229", "002", "003-006" 10 20;
%                   ]

%% load data
% first load from ret
% extract goodfit cells, keep track of total #
% store variable for stim to RF center distance
% then load size tuning data from only good fit cells
% think total distance OK
% stimulus is vertical, could consider x and y separately

nrun = size(dataSelect,1)
% load retinotopy data
fprintf('\nLoading retinotopy data:\n')
% need just the good fit indices, RF centers
retData = cell(nrun, 3); % 3 columns for goodfit_ind, cellAz, cellEl
for irun = 1:nrun
    date = char(dataSelect(irun,2));
    mouse = char(dataSelect(irun,1));
    ret_str = char(dataSelect(irun,3));
    RetFolder = [date '_' mouse '_runs-' ret_str];
    fprintf(['irun: ' num2str(irun) '/' num2str(nrun) ' - ret folder: ' RetFolder '\n'])
    
    %CD = ['/Users/kevinmurgas/Documents/MATLAB/Analysis/Analysis/2P/' date '_' mouse '/' RetFolder];
    CD = ['H:\home\kevin\Analysis\2P\' date '_' mouse '\' RetFolder];
    %CD = ['H:\home\valerie\Data\2p\' date '_' mouse '\' RetFolder];
    %CD = ['H:\home\lindsey\Data\2P_images\' date '_' mouse '\' RetFolder];
    cd(CD);
    
    % load ret data:
    % specifically RF centers
    % []_lbub_fits.mat contains: goodfit_ind, lbub_diff, lbub_fits, resp_ind
    % goodfit_ind describes cells to extract
    % lbub_fits contains RF fit data
    fprintf('Loading retinotopy fit data\n')
    fn_out = [RetFolder '_lbub_fits.mat'];
    load(fn_out);
    retData{irun,1} = goodfit_ind;
    retData{irun,2} = lbub_fits(goodfit_ind,4,4); % cellAz
    retData{irun,3} = lbub_fits(goodfit_ind,5,4); % cellEl
end
fprintf('Retinotopy fits loaded\n')

% load size tuning data
fprintf('\nLoading size tuning data:\n')
szData = cell(nrun, 4);
% 4 columns for sizeTune (nSize, nCon, nCells), sizeMean, sizeSEM, cellDists
for irun = 1:nrun
    date = char(dataSelect(irun,2));
    mouse = char(dataSelect(irun,1));
    sz_str = char(dataSelect(irun,4));
    SzFolder = [date '_' mouse '_runs-' sz_str];
    fprintf(['irun: ' num2str(irun) '/' num2str(nrun) ' - ret folder: ' SzFolder '\n'])
    
    %CD = ['/Users/kevinmurgas/Documents/MATLAB/Analysis/Analysis/2P/' date '_' mouse '/' SzFolder];
    CD = ['H:\home\kevin\Analysis\2P\' date '_' mouse '\' SzFolder];
    %CD = ['H:\home\valerie\Data\2p\' date '_' mouse '\' SzFolder];
    %CD = ['H:\home\lindsey\Data\2P_images\' date '_' mouse '\' SzFolder];
    cd(CD);
    
    % load size tuning data:
    % first load some input data, variables
    load([SzFolder '_input.mat'])
    load([SzFolder '_Tuning.mat'],'szs')
    
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
    
    % specifically sizeTune and cellDists from sizeTuneData.mat
    % add catch for if sizeTuneData.mat does not exist, to then compute
    fprintf('Loading size tuning data\n')
    fn_out = [SzFolder '_sizeTuneData.mat'];
    if exist(fn_out, 'file')
        load(fn_out);
        szData{irun,1} = sizeTune;
        szData{irun,2} = sizeMean;
        szData{irun,3} = sizeSEM;
        szData{irun,4} = cellDists;
    else
        fprintf('sizeTuneData.mat not found, computing and saving\n')
        
        load([SzFolder '_Tuning.mat'])
        
        goodfit_ind = retData{irun,1};
        stimEl = double(input.gratingElevationDeg);
        stimAz = double(input.gratingAzimuthDeg);
        cellDists = sqrt((retData{irun,2}-stimAz).^2+(retData{irun,3}-stimEl).^2);
        
        nCells = length(goodfit_ind);
        sizeTune = cell(nSize,nCon,nCells);
        sizeMean = zeros(nSize,nCon,nCells);
        sizeSEM = sizeMean;
        
        for i = 1:length(goodfit_ind)
            iCell = goodfit_ind(i);
            
            for iCon = 1:nCon
                for iSize = 1:nSize
                    ind_all = intersect(Ind_struct(iSize).all_trials,conInds{iCon});
                    %stimOff = mean(mean(tc_dfof((nOff/2):nOff,iCell,ind_all),3),1);
                    
                    % this takes mean dF/F during stimOn at given ind_all
                    stimOn = mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),1);
                    sd = std(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),1));
                    sizeMean(iSize,iCon,i) = mean(stimOn,3);
                    sizeSEM(iSize,iCon,i) = sd/sqrt(length(ind_all));
                    
                    % this cell matrix stores all stimOn for all trials in ind_all
                    % dims (szs, cons, cells)
                    sizeTune{iSize, iCon, i} = squeeze(stimOn);
                end
            end
        end
        % save
        filename = [SzFolder '_sizeTuneData.mat'];
        save(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
        fprintf('Saved sizeTuneData.mat\n')
        
        szData{irun,1} = sizeTune;
        szData{irun,2} = sizeMean;
        szData{irun,3} = sizeSEM;
        szData{irun,4} = cellDists;
    end
end
fprintf('Size tuning data loaded\n')

% load fit data
fprintf('\nLoading size tuning fits:\n')
fitData = cell(nrun, 3);
% 3 columns for allFtest, allPrefSize, allSuppInd
for irun = 1:nrun
    date = char(dataSelect(irun,2));
    mouse = char(dataSelect(irun,1));
    sz_str = char(dataSelect(irun,4));
    SzFolder = [date '_' mouse '_runs-' sz_str];
    fprintf(['irun: ' num2str(irun) '/' num2str(nrun) ' - ret folder: ' SzFolder '\n'])
    
    %CD = ['/Users/kevinmurgas/Documents/MATLAB/Analysis/Analysis/2P/' date '_' mouse '/' SzFolder];
    CD = ['H:\home\kevin\Analysis\2P\' date '_' mouse '\' SzFolder];
    %CD = ['H:\home\valerie\Data\2p\' date '_' mouse '\' SzFolder];
    %CD = ['H:\home\lindsey\Data\2P_images\' date '_' mouse '\' SzFolder];
    cd(CD);
    
    % load fit data:
    %load([SzFolder '_sizeFitResults.mat'])
    load([SzFolder '_sizeFitResults_SP.mat'])
    % loads sizeFits struct containing:
    % nCellsxnCon - Rsq1, Rsq2, prefSize, suppInd, Fscore, Ftest
    % nCellsxnConx3 - coefs1
    % nCellsxnConx6 - coefs2
    
    fitData{irun,1} = sizeFits.Ftest;
    fitData{irun,2} = sizeFits.prefSize;
    fitData{irun,3} = sizeFits.suppInd;
end
fprintf('Size fitting data loaded\n')

% concatenate all runs?


%% present data
% ideas:
% mean size tuning curve for each cell? prob too many, select certain cells
% mean size tuning curve for each contrast
% histograms of stim to RF center distance for bins
% mean size tuning curve by bins
% rough estimates of pref size (peak), SI ([peak-largest]/peak)
% ^ by contrast and by bins
%
% after adding fit data:
% fit prefSize and SI, histograms and by contrast
% by contrast overall, by distance bins
%
% look at valeries data

% plot mean and raw size tuning for a selected cell and con
% figure(1)
% iCell = 6;
% iCon = 2;
% errorbar(szs,sizeMean(:,iCon,iCell),sizeSEM(:,iCon,iCell))
% hold on
% sizeX = [];
% sizeY = [];
% for iSz = 1:length(szs)
%     sizeX = [sizeX szs(iSz).*ones(1,length(sizeTune{iSz,iCon,iCell}))];
%     sizeY = [sizeY sizeTune{iSz,iCon,iCell}'];
% end
% plot(sizeX, sizeY, '.b')

% first concatenate all sizeMean, sizeSEM along 3rd dimension
fprintf('\nConcatenating size tuning data across runs\n')
catSzMean = [];
catSzSEM = [];
catCellDists = [];
catFtest = [];
catPrefSize = [];
catSuppInd = [];
Syt7KOflag = [];
for irun = 1:nrun
    catSzMean = cat(3,catSzMean,szData{irun,2});
    catSzSEM = cat(3,catSzSEM,szData{irun,3});
    catCellDists = cat(1,catCellDists,szData{irun,4});
    
    catFtest = cat(1,catFtest,fitData{irun,1});
    catPrefSize = cat(1,catPrefSize,fitData{irun,2});
    catSuppInd = cat(1,catSuppInd,fitData{irun,3});
    
    if dataSelect(irun,1)=="i927"
        Syt7KOflag = cat(1,Syt7KOflag,ones(size(fitData{irun,1},1),1));
    else
        Syt7KOflag = cat(1,Syt7KOflag,zeros(size(fitData{irun,1},1),1));
    end
    
    fprintf(['irun: ' num2str(irun) '/' num2str(nrun) ', with ' num2str(length(szData{irun,4})) ' cells\n'])
end
nCells = length(catCellDists);
Syt7KOflag = logical(Syt7KOflag);
nWT = sum(~Syt7KOflag);
nKO = sum(Syt7KOflag);
suppCells = find(catFtest);

% for control
catSuppInd(find(catSuppInd<0))=0;
catSuppInd(find(catSuppInd>1))=1;

figure(2);clf;
histogram(catCellDists,30)
title(['RF-stim distances, n=' num2str(nCells)])
xlabel('Distance (deg)')

%% tuning curve stuff
% *** separate supp vs non supp
% mean size tuning curve for 4 contrasts, average all cells
figure(3);clf;
% one plot, 4 lines
subplot(2,1,1)
for iCon = 1:nCon
    errorbar(szs,mean(catSzMean(:,iCon,~Syt7KOflag),3),1/sqrt(nWT)*geomean(catSzSEM(:,iCon,~Syt7KOflag),3))
    hold on
end
title(['Average Size Tuning Curve WT, n=' num2str(nWT)])
xlabel('Size (deg)')
ylabel('dF/F')
legend(num2str(cons'))
subplot(2,1,2)
for iCon = 1:nCon
    errorbar(szs,mean(catSzMean(:,iCon,Syt7KOflag),3),1/sqrt(nKO)*geomean(catSzSEM(:,iCon,Syt7KOflag),3))
    hold on
end
title(['Average Size Tuning Curve Syt7 KO, n=' num2str(nKO)])
xlabel('Size (deg)')
ylabel('dF/F')
legend(num2str(cons'))

% do binning, based on cellDists
fprintf('\nSorting cells by RF center distance\n')
bIndex = 0*catCellDists;
nBin = 3;
for i=1:nCells
    if catCellDists(i)<6
        bIndex(i) = 1;
    elseif catCellDists(i)<12
        bIndex(i) = 2;
    else
        bIndex(i) = 3;
    end
end
nBinCells = [sum(bIndex==1) sum(bIndex==2) sum(bIndex==3)];
fprintf([num2str(nBinCells) ' cells in each bin\n'])

binStrs = ["0-6","6-12","12+"];
cBins = categorical({'0-6','6-12','12+'});
cBins = reordercats(cBins,{'0-6','6-12','12+'});
fprintf(char(binStrs))
fprintf('\n')

% cellDists histogram by bin
figure(4);clf;
bar(cBins,nBinCells)
title('Cell Counts by RF Distance Bin')
xlabel('Distance Bin')
ylabel('# cells')

% *** add WT vs Syt7, non-supp vs supp ***
% mean size tuning curve for 4 contrasts, separate by bins
figure(5);clf;
for k = 1:nBin
    subplot(1,nBin,k)
    for iCon=1:nCon
        errorbar(szs,mean(catSzMean(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(catSzSEM(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['Average Size Tuning Curve, RF dist ' char(binStrs(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    %ylim([0 0.08])
    legend(num2str(cons'))
end

cCons = categorical(cons);

% **** expand this to show plots separate for WT and Syt7 ***
% add a plot for the number of unsuppressed vs suppressed cells
figure(6);clf;
for k=1:nBin
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    subplot(1,3,k)
    
    binFtest = zeros(nCon,2); % one column each for unsupp, supp cells
    if nCutCells
        for iCon = 1:nCon
            binFtest(iCon,1) = nCutCells - sum(catFtest(binCells,iCon));
            binFtest(iCon,2) = sum(catFtest(binCells,iCon));
        end
    end
    bar(cCons,binFtest)
    title(['Nonsupp vs Supp cells, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('# cells')
    legend('Nonsupp','Supp')
end

% *** change this, make it show all three bins ***
nCutCells = zeros(nCon,1);
% plot pref size and SI for WT vs Syt7
figure(7);clf;
wtCells = find((catCellDists<15).*(~Syt7KOflag));
% first show average PrefSize at each contrast only suppressed cells
for iCon = 1:nCon
    suppCellInds = find(catFtest(wtCells,iCon));
    nCutCells(iCon) = length(suppCellInds);
    prefAvWT(iCon) = mean(catPrefSize(suppCellInds,iCon),1);
    prefSEMWT(iCon) = std(catPrefSize(suppCellInds,iCon))/sqrt(nCutCells(iCon));
    SIavWT(iCon) = mean(catSuppInd(suppCellInds,iCon),1);
    SISEMWT(iCon) = std(catSuppInd(suppCellInds,iCon))/sqrt(nCutCells(iCon));
end
koCells = find((catCellDists<15).*(Syt7KOflag));
% first show average PrefSize at each contrast only suppressed cells
for iCon = 1:nCon
    suppCellInds = find(catFtest(koCells,iCon));
    nCutCells(iCon) = length(suppCellInds);
    prefAvKO(iCon) = mean(catPrefSize(suppCellInds,iCon),1);
    prefSEMKO(iCon) = std(catPrefSize(suppCellInds,iCon))/sqrt(nCutCells(iCon));
    SIavKO(iCon) = mean(catSuppInd(suppCellInds,iCon),1);
    SISEMKO(iCon) = std(catSuppInd(suppCellInds,iCon))/sqrt(nCutCells(iCon));
end
subplot(1,2,1)
bar(cCons,[prefAvWT;prefAvKO]','grouped')
% hold on
% errorbar(cCons,prefAvWT,prefSEMWT,'.')
% hold off
title('Preferred Size, RF dist<15')
xlabel('Contrast')
ylabel('Average Pref Size')
ylim([0 50])
legend(['WT(n=' num2str(length(wtCells)) ' cells)'],['Syt7 KO (n=' num2str(length(koCells)) ' cells)'])
% second show average SI at each contrast
subplot(1,2,2)
bar(cCons,[SIavWT;SIavKO]','grouped')
% hold on
% errorbar(cCons,SIavWT,SISEMWT,'.')
% hold off
title('Suppression Index, RF dist<15')
xlabel('Contrast')
ylabel('Average Supp Ind')
ylim([0 0.8])
legend('WT','Syt7 KO')

% % rough estimates, scrapped
% % calculate rough estimates of prefSize and SI
% % go through each cell, measure rough prefsize and SI at each con
% roughPrefSize = zeros(nCells, nCon);
% roughSI = zeros(nCells, nCon);
% for iCell = 1:nCells
%     for iCon = 1:nCon
%         prefInd = find(catSzMean(:,iCon,iCell) == max(catSzMean(:,iCon,iCell)));
%         roughPrefSize(iCell,iCon) = szs(prefInd);
%         roughSI(iCell,iCon) = 1 - catSzMean(end,iCon,iCell)/catSzMean(prefInd,iCon,iCell);
%     end
% end
% 
% % present prefSize and SI by contrast, separate by bins
% 
% figure(16);clf;
% for k=1:3
%     goodCutCells = find(bIndex==k);
%     nCutCells = length(goodCutCells);
%     
%     % first show average PrefSize at each contrast
%     subplot(2,3,k)
%     prefAv = mean(roughPrefSize(goodCutCells,:),1);
%     prefSEM = std(roughPrefSize(goodCutCells,:),0,1)/sqrt(nCutCells);
%     bar(cCons,prefAv)
%     hold on
%     errorbar(cCons,prefAv,prefSEM,'.')
%     hold off
%     title(['Rough Preferred Size, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
%     xlabel('Contrast')
%     ylabel('Average Pref Size')
%     ylim([0 60])
%     
%     % second show average SI at each contrast
%     subplot(2,3,k+3)
%     SIav = mean(roughSI(goodCutCells,:),1);
%     SISEM = std(roughSI(goodCutCells,:),0,1)/sqrt(nCutCells);
%     bar(cCons,SIav)
%     hold on
%     errorbar(cCons,SIav,SISEM,'.')
%     hold off
%     title(['Rough Suppression Index, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
%     xlabel('Contrast')
%     ylabel('Average Supp Ind')
%     ylim([0 0.6])
% end

% cumulative distributions of prefSize and SI, with contrasts overlaid
%prefSize
figure(8);clf;
for k=1:nBin
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    binWTcells = intersect(binCells,find(~Syt7KOflag));
    nWTcells = length(binWTcells);
    binKOcells = intersect(binCells,find(Syt7KOflag));
    nKOcells = length(binKOcells);
    
    % row 1=WT
    subplot(2,3,k)
    for iCon=1:nCon
        %cdfplot(catPrefSize(binWTcells,iCon)) %both supp and nonsupp cells
        binWTcells_supp = intersect(binWTcells,find(catFtest(:,iCon)));
        cdfplot(catPrefSize(binWTcells_supp,iCon))
        hold on
        xlim([0 100])
    end
    hold off
    title(['WT, RF dist: ' char(binStrs(k)) ' (n=' num2str(nWTcells) ' cells)'])
    xlabel('Pref Size')
    ylabel('frac cells')
    legend(num2str(cons'),'location','se')
    % row 2=KO
    subplot(2,3,k+3)
    for iCon=1:nCon
        %cdfplot(catPrefSize(binKOcells,iCon)) %both supp and nonsupp cells
        binKOcells_supp = intersect(binKOcells,find(catFtest(:,iCon)));
        cdfplot(catPrefSize(binKOcells_supp,iCon))
        hold on
        xlim([0 100])
    end
    hold off
    title(['Syt7KO, RF dist: ' char(binStrs(k)) ' (n=' num2str(nKOcells) ' cells)'])
    xlabel('Pref Size')
    ylabel('frac cells')
    legend(num2str(cons'),'location','se')
end

%SI
figure(9);clf;
for k=1:nBin
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    binWTcells = intersect(binCells,find(~Syt7KOflag));
    nWTcells = length(binWTcells);
    binKOcells = intersect(binCells,find(Syt7KOflag));
    nKOcells = length(binKOcells);
    
    % row 1=WT
    subplot(2,3,k)
    for iCon=1:nCon
        %cdfplot(catSuppInd(binWTcells,iCon)) %both supp and nonsupp cells
        binWTcells_supp = intersect(binWTcells,find(catFtest(:,iCon)));
        cdfplot(catSuppInd(binWTcells_supp,iCon))
        hold on
        xlim([-1 2])
    end
    hold off
    title(['WT, RF dist: ' char(binStrs(k)) ' (n=' num2str(nWTcells) ' cells)'])
    xlabel('SI')
    ylabel('frac cells')
    legend(num2str(cons'),'location','nw')
    % row 2=Syt7KO
    subplot(2,3,k+3)
    for iCon=1:nCon
        %cdfplot(catSuppInd(binKOcells,iCon)) %both supp and nonsupp cells
        binKOcells_supp = intersect(binKOcells,find(catFtest(:,iCon)));
        cdfplot(catSuppInd(binKOcells_supp,iCon))
        hold on
        xlim([-1 2])
    end
    hold off
    title(['Syt7KO, RF dist: ' char(binStrs(k)) ' (n=' num2str(nKOcells) ' cells)'])
    xlabel('SI')
    ylabel('frac cells')
    legend(num2str(cons'),'location','nw')
end

%% Fit prefSize changes with contrast
% will need to store fit params for examining a single size

% normalize by 0.8 con pref size
normPrefSize = catPrefSize./catPrefSize(:,nCon);

% first calculate slopes, fitting with a 1st degree polynomial (line)
lineFits = zeros(nCells,2);
for i=1:nCells
    lineFits(i,:) = polyfit(cons,normPrefSize(i,:),1);
end

% the plot each cell's prefSize vs con curve, within each bin
for k=1:nBin
    figure(8+k);clf;
    [n, n2] = subplotn(nBinCells(k));
    binInds = find(bIndex == k);
    for i = 1:nBinCells(k)
        iCell = binInds(i);
        subplot(n,n2,i)
        plot(cons,normPrefSize(iCell,:),'o')
        hold on
        plot(cons,polyval(lineFits(iCell,:),cons),'r-')
        hold off
        title(['Cell #' num2str(iCell) ', m=' num2str(lineFits(iCell,1))])
        xlabel('Contrast')
        ylabel('Pref Sz')
        ylim([0 2])
    end
end

% readout means and stdevs for each bin
for k=1:nBin
    fprintf(['Bin #' num2str(k) ', Mean slopes: ' num2str(mean(lineFits(bIndex==k,1),1)) ', stdev: ' num2str(std(lineFits(bIndex==k,1),[],1)) '\n'])
end
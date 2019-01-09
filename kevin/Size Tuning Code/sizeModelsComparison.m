%% this script was used to compared two models for size tuning fits
% full model uses sigmoids fitting with centers (6 fit param)
% reduced model sets centers to zero (4 fit param)
% comparing by sequential F-test on each condition (cell x Con)
% also uses Akaike's Information Criterion (AIC)
%
% no longer using script... decided on reduced model

%% Load data
% currently need to run sizeTuningAfterRet, runs through only one dataset
%
% insert here code to load/define variables including:
% good fit indices
% fit data at good fits (RF center), stim Az+El, select cells within range
% size + con ranges (szs)? need to define size tuning curve
% size tuning data
% can copy some of of this code over
loadFlag = 0;
if loadFlag
    % load tc data
    
    %load ret data
    % load ret fits - loads 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind'
    fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
    fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
    load(fn_out);
    cellAz = lbub_fits(:,4,4);
    cellEl = lbub_fits(:,5,4);
    fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
end

% input stimulus location based on experimental choice
stimEl = -15;
stimAz = 5;
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])


fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz-stimAz).^2+(cellEl-stimEl).^2);

% compare distance to select cells
cutOffRadius = 7;
fprintf(['Isolating cells with RF centers within ' num2str(cutOffRadius) ' degrees\n'])

cutCells = find(cellDists < 7)';
goodCutCells = intersect(cutCells,goodfit_ind);
% goodCutCells = setdiff(cutCells,goodfit_ind); % set difference, actually bad fit cells
% cutCells = [49 111];
nCutCells = length(goodCutCells);
fprintf([num2str(nCutCells) ' cells selected\n'])


conTrials = cell2mat(input.tGratingContrast);
cons = unique(conTrials);
nCon = length(cons);
conInds = cell(nCon,1);
for i = 1:nCon
    conInds{i} = find(conTrials == cons(i));
end

sizeTune = zeros(nStim,nCon,nCutCells);
sizeSEM = sizeTune;
[n, n2] = subplotn(nCutCells);
figure(3);clf;
suptitle(['Size Tuning Curves for cells within ' num2str(cutOffRadius) ' deg of stim'])
for i = 1:nCutCells
    iCell = goodCutCells(i);
    % this first section shows average trials for each stim
    % should show average of each of three trials
    %     figure;
    %     for iCond = 1:nStim
    %         subplot(4,2,iCond)
    %         ind_all = Ind_struct(iCond).all_trials;
    %         plot(squeeze(mean(tc_dfof(:,iCell,ind_all),3)))
    %         title(['Size: ' num2str(szs(iCond)) ' deg'])
    %         ylim([-0.05 0.4])
    %     end
    %     suptitle(['Cell #: ' num2str(iCell)])
    
    % this section shows size tuning curves
    subplot(n,n2,i)
    for iCon = 1:nCon
        for iSize = 1:nStim
            ind_all = intersect(Ind_struct(iSize).all_trials,conInds{iCon});
            %stimOff = mean(mean(tc_dfof((nOff/2):nOff,iCell,ind_all),3),1);
            stimOn = mean(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),3),1);
            sd = std(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind_all),1));
            sizeTune(iSize,iCon,i) = stimOn;
            sizeSEM(iSize,iCon,i) = sd/sqrt(length(ind_all));
        end
        errorbar(szs,sizeTune(:,iCon,i),sizeSEM(:,iCon,i))
        hold on
    end
    ylim([-1.2*(max(max(sizeSEM(:,:,i)))) 1.2*(max(max(sizeTune(:,:,i)))+max(max(sizeSEM(:,:,i))))])
    title(['Cell #: ' num2str(iCell)])
    hold off
end
% Construct a Legend with the data from the sub-plots
hL = legend(num2str(cons'));
% Programatically move the Legend
newPosition = [0.5 0.07 0.2 0.2];
newUnits = 'normalized';
set(hL, 'Position', newPosition, 'Units', newUnits);


%% Fit size tuning curves
chosen = 94;

% full model, with parameters for sigmoid centers
allRsq1 = zeros(nCutCells, nCon); % cells x contrasts, collects rsq for each fit
allSSE1 = allRsq1; % also collecting SSE for each fit

allPrefSize1 = allRsq1; % these two just collect parameters of interest
allSuppInd1 = allRsq1;

figure(1);clf;
for i = 1:nCutCells
    iCell = goodCutCells(i);
    fprintf(['Cell# ' num2str(iCell)])
    dummyTune = sizeTune(:,:,i);
    dummySEM = sizeSEM(:,:,i);
    szRng = linspace(0,max(szs));
    
    %logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
    logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-coefs(6))))
    logfit3 = @(coefs,xdata) coefs(1)./(1+exp(-1.38*(xdata-25)))
    
    for iCon = 1:nCon
        dumdum = [0 dummyTune(:,iCon)'];
        szs0 = [0 szs];
        SStot = sum((dumdum-mean(dumdum)).^2);
        
        guess = [2*max(dumdum) 0.35 5 mean(dumdum) 0.06 10]; % first curve should be steeper
        lb = [0 0 0 0 0 0];
        ub = [inf inf max(szs) inf inf 2*max(szs)];
        [fitx,resnorm] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub);
        r2 = 1 - resnorm/SStot; % resnorm is SSE
        allRsq1(i, iCon) = r2;
        allSSE1(i, iCon) = resnorm;
        
        fprintf(['Fit R-sq:' num2str(r2) '\n\n'])
        
        fprintf('Selecting fit 2\n')
        fitout = logfit2(fitx,szRng);
        
        fprintf('Fit readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ' centered at: ' num2str(fitx(3))])
        fprintf(['\nInh amp: ' num2str(fitx(4)) ' centered at: ' num2str(fitx(6))])
        
        maxResp = max(fitout);
        prefSize = szRng(find(fitout==maxResp,1));
        fprintf(['\nPref size: ' num2str(prefSize) ' (dF/F: ' num2str(maxResp) ')\n'])
        allPrefSize1(i, iCon) = prefSize;
        
        suppInd = 1 - fitout(end)/maxResp;
        fprintf(['Suppression Index: ' num2str(suppInd) '\n'])
        allSuppInd1(i, iCon) = suppInd;
        
        if iCell == chosen % only plot first cell
            subplot(1,nCon,iCon)
            errorbar(szs0,dumdum,[0 dummySEM(:,iCon)'])
            hold on
            plot(szRng,fitout,'-')
            hold off
            ylim([0 1.4*maxResp])
            title(['Cell #' num2str(iCell) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(r2) ', SI=' num2str(suppInd)])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
        end
    end
end

% reduced model, with zero center
allRsq2 = zeros(nCutCells, nCon); % cells x contrasts, collects rsq for each fit
allSSE2 = allRsq2; % also collecting SSE for each fit

allPrefSize2 = allRsq2; % these two just collect parameters of interest
allSuppInd2 = allRsq2;

figure(2);clf;
for i = 1:nCutCells
    iCell = goodCutCells(i);
    fprintf(['Cell# ' num2str(iCell)])
    dummyTune = sizeTune(:,:,i);
    dummySEM = sizeSEM(:,:,i);
    szRng = linspace(0,max(szs));
    
    logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata))) - coefs(3)./(1+exp(-coefs(4)*(xdata)))
    
    for iCon = 1:nCon
        dumdum = [0 dummyTune(:,iCon)'];
        szs0 = [0 szs];
        SStot = sum((dumdum-mean(dumdum)).^2);
        
        guess = [2*max(dumdum) 0.35 mean(dumdum) 0.06]; % first curve should be steeper
        lb = [0 0 0 0];
        ub = [inf inf inf inf];
        [fitx,resnorm] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub);
        r2 = 1 - resnorm/SStot; % resnorm is SSE
        allRsq2(i, iCon) = r2;
        allSSE2(i, iCon) = resnorm;
        
        fprintf(['Fit R-sq: ' num2str(r2) '\n\n'])
        
        fprintf('Selecting fit 2\n')
        fitout = logfit2(fitx,szRng);
        
        fprintf('Fit readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ' centered at 0'])
        fprintf(['\nInh amp: ' num2str(fitx(3)) ' centered at 0'])
        
        maxResp = max(fitout);
        prefSize = szRng(find(fitout==maxResp,1));
        fprintf(['\nPref size: ' num2str(prefSize) ' (dF/F: ' num2str(maxResp) ')\n'])
        allPrefSize2(i, iCon) = prefSize;
        
        suppInd = 1 - fitout(end)/maxResp;
        fprintf(['Suppression Index: ' num2str(suppInd) '\n'])
        allSuppInd2(i, iCon) = suppInd;
        
        if iCell == chosen % only plot first cell
            subplot(1,nCon,iCon)
            errorbar(szs0,dumdum,[0 dummySEM(:,iCon)'])
            hold on
            plot(szRng,fitout,'-')
            hold off
            ylim([0 1.4*max([maxResp max(dumdum)])])
            title(['Cell #' num2str(iCell) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(r2) ', SI=' num2str(suppInd)])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('data','fit')
        end
    end
end

%% examine parameters
figure(3);clf;
subplot(1,3,1)
histogram(allPrefSize1(:),50)
ylim([0 60])
title('PrefSize for full model')
subplot(1,3,2)
histogram(allPrefSize2(:),50)
ylim([0 60])
title('PrefSize for reduced model')
subplot(1,3,3)
histogram(allPrefSize1(:)-allPrefSize2(:),50)
ylim([0 60])
title('Paired Difference (full-reduced)')

figure(4);clf;
subplot(1,3,1)
histogram(allSuppInd1(:),30)
title('SuppInd for full model')
ylim([0 60])
subplot(1,3,2)
histogram(allSuppInd2(:),30)
ylim([0 60])
title('SuppInd for reduced model')
subplot(1,3,3)
histogram(allSuppInd1(:)-allSuppInd2(:),30)
ylim([0 60])
title('Paired Difference (full-reduced)')

% allSSE1 refers to full model, allSSE2 refers to reduced model
figure(5);clf;
subplot(1,3,1)
histogram(allRsq1(:),30)
ylim([0 20])
title('R^2 for full model')
subplot(1,3,2)
histogram(allRsq2(:),30)
ylim([0 20])
title('R^2 for reduced model')
subplot(1,3,3)
histogram(allRsq1(:)-allRsq2(:),30)
ylim([0 20])
title('Paired Difference (full-reduced)')

% insert adjusted Rsq *********

%% Do F-test and calculate AIC
% F-test
dfR = 9 - 4; % reduced model df = #sizes - 4 model parameters
dfF = 9 - 6; % full model df = #sizes - 6 model parameters
allFscore = ((allSSE2(:) - allSSE1(:))/(dfR-dfF)) ./ (allSSE1(:)/dfF);
Fcrit = 5.4095; % cutoff for df=3 and df=5
% Fcrit = 9.0135; % cutoff for df=5 and df=3
portionReject = sum(allFscore > Fcrit)/length(allFscore)

figure(6);clf;
histogram(allFscore,50)
hold on
line([Fcrit Fcrit],[0 100],'Color','red')
text(20,40,['Portion reject (\alpha = 0.05): ' num2str(portionReject)])
hold off
title('Sequential F-score for Full Model vs Reduced Model')

% AIC
% N*ln(SS/N) + 2K where N= number of data points (9), K = parameters + 1
% correction term: + 2K(K+1)/(N-K-1)
allAICa = 9*log(allSSE2/9) + 10; % reduced model K=5
allAICb = 9*log(allSSE1/9) + 14; % reduced model K=7
AICCa = allAICa + 2*5*6/3; %+20
AICCb = allAICb + 2*7*8; %+112
delAIC = allAICb - allAICa;
delAICc = AICCb-AICCa;
% probability = exp(-0.5del)/(1+exp(-0.5del))
% these are probabilities that model B (full, SSE1) is true
prob = exp(-0.5*delAIC)./(1+exp(-0.5*delAIC));
fracFull = sum(prob(:)>0.5)/length(prob(:))
probC = exp(-0.5*delAICc)./(1+exp(-0.5*delAICc));
fracFullC = sum(probC(:)>0.5)/length(probC(:))

figure(7);clf;
histogram(delAIC,30)
hold on
text(-30,14,['Portion >0 (full model prevails): ' num2str(fracFull)])
hold off
title('\Delta AIC for Full Model - Reduced Model')

figure(8);clf;
histogram(delAICc,30)
hold on
text(65,14,['Portion >0 (full model prevails): ' num2str(fracFullC)])
hold off
title('\Delta AICc for Full Model - Reduced Model')

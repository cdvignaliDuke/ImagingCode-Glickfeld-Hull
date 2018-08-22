close all
clear all

FSAV_V1_decode

load('Z:\Analysis\FSAV Choice\FSAV_decodeData.mat')
load('Z:\Analysis\FSAV Summaries\FSAV_V1_decode\V1_decode_startAlign_FSAV_attnData')
fnout = 'Z:\Analysis\FSAV Choice';
bxParams_FSAV
nexp = length(dcExpt);

doExptPlots = 0;
%%
outcomecolmat = {'k';'r'};
avcolmat = {'k','c'};

ampBinsFine = [0 exp(linspace(log(0.0001),log(1),7))];
ampBins = [0 0.0001 0.1 1];
weightsDiffZeroAlpha = 0.05/2;
nPC = 15;
weightLim = [-3 4];
weightSubLim = [-2.5 2.5];
corrLim = [-1 1];
siLim = [-12 12];
adaptLim = [-1.5 1.5];
maxCellN = 10;
nBoot = 1000;

minTrN = 20;
theta90Threshold = 11.25;
decisionVariable = 0.5;
detectThresholdHR = 0.8;
detectThresholdOri = 45;

oriBinsFine = [0 1 16 23 32 45 64 90];
oriBins = [0 1 32 90];
%%

amplitudes = [];
orientations = [];
invOrientations = [];
for iexp = 1:nexp
    tOri = dcExpt(iexp).trialOrientation;
    tOri = binAndRelabelTrialOrientation(tOri,oriBinsFine);
    orientations = cat(2,orientations,unique(tOri));
    tAmp = dcExpt(iexp).audTrialAmp;
    tAmp = binAndRelabelTrialOrientation(tAmp,ampBinsFine);
    amplitudes = cat(2,amplitudes,unique(tAmp));
    if ~isempty(dcExpt(iexp).catchResp)
        tOri = binAndRelabelTrialOrientation(...
            dcExpt(iexp).catchOrientation,oriBinsFine);
        invOrientations = cat(2,invOrientations,unique(tOri));
    end
end
allOris = unique(orientations);
allAmps = unique(amplitudes);
allInvOris = unique(invOrientations);

%%
targetRespFig = figure;
respCells_target = cell(1,nexp);
respCells_base = cell(1,nexp);
valRespOri = nan(3,4,nexp);
invRespOri = nan(3,4,nexp);
catchInd = false(1,nexp);
for iexp = 1:nexp
if ~isempty(dcExpt(iexp).catchResp)        
        tInvOri = binAndRelabelTrialOrientation(...
            dcExpt(iexp).catchOrientation,oriBins);
        invOrientations = unique(tInvOri);     
%         invRespOriExp
        if any(invOrientations ~= 0)  

            
            respCells_target{iexp} = logical(...
                dcExpt(iexp).targetResponsiveCells);
            respCells_base{iexp} = logical(...
                dcExpt(iexp).firstBaseResponsiveCells);
            
            valResp = dcExpt(iexp).stimResp';
            invResp = dcExpt(iexp).catchResp';
            tOri = binAndRelabelTrialOrientation(...
                dcExpt(iexp).trialOrientation,oriBins);
            nori = length(invOrientations);
            ncells = size(valResp,2);
            
            targetOnly = respCells_target{iexp} & ~respCells_base{iexp};
            targetAndBase = respCells_target{iexp} & respCells_base{iexp};
            baseOnly = ~respCells_target{iexp} & respCells_base{iexp};

            cellTypeInd = {targetOnly,targetAndBase,baseOnly,...
                targetOnly | targetAndBase | baseOnly};
            cellTypeName = {'Target Only';'Target & Base';'Base Only';'All'};

            invRespMean = nan(nori,4);
            valRespMean = nan(nori,4);
            for iori = 1:nori
                ind1 = tInvOri == invOrientations(iori);
                n = sum(ind1);
                ind2 = randsample(find(tOri == invOrientations(iori)),n);
                for i = 1:4
                    cellInd = cellTypeInd{i};
                    invRespMean(iori,i) = mean(mean(invResp(ind1,cellInd),2),1);
                    valRespMean(iori,i) = mean(mean(valResp(ind2,cellInd),2),1);
                end
            end

            figure(targetRespFig)
%             suptitle([dcExpt(iexp).mouse '-' dcExpt(iexp).date])
            x = invOrientations;
            for i = 1:4
                subplot(1,4,i)
                hold on
                y = valRespMean(:,i) - invRespMean(:,i);
                h = plot(x,y,'-');
                h.Color = cueColor{valid};
%                 y = invRespMean(:,i);
%                 h = plot(x,y,'-');
%                 h.Color = cueColor{invalid};
                title(cellTypeName{i})
                figXAxis([],'Target Orientation (deg)',[-10 100],invOrientations,{'0';'16-32';'32-90'});
                figYAxis([],'Valid - Invalid dF/F (%)',[-0.07 0.12])
                figAxForm
                valRespOri(:,i,iexp) = valRespMean(:,i);
                invRespOri(:,i,iexp) = invRespMean(:,i);
            end
           
            catchInd(iexp) = true;
        end
end
end
figure(targetRespFig)
for i = 1:4
    subplot(1,4,i)
    y = nanmean(valRespOri(:,i,catchInd) - invRespOri(:,i,catchInd),3);
    yerr = ste(valRespOri(:,i,catchInd) - invRespOri(:,i,catchInd),3);
    h = errorbar(x,y,yerr,'ko-');
    h.LineWidth = 2;   
    hline(0,'k--')
end

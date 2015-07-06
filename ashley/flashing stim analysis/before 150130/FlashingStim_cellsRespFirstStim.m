mouse = 'AW07';
SubNum = '607';
date = '141215';
ImgFolder = '004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructVar.mat');
load('dataTC.mat');
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% dirTuning dataset
%get average response to vis stim for each cycle
siz = size(dataStructVar.Cycles,2);
data = dataTC.dirTuning;
nCells = size(data{1},2);

for icyc = 1:siz
    for icell = 1:nCells
        thisL = dataStructVar.cycTrialL{icyc};
        thisData = data{icyc};
        for itrial = 1:dataStructVar.cycNTrials{icyc}
            thisMat(:,itrial,icell) = thisData(1+(thisL.*(itrial-1)):thisL.*itrial,icell);
        end
        dataTrialsCells{icyc} = thisMat;
    end
    clear thisL thisData thisMat
end

start = 1;
for ilength = dataStructVar.Cycles
    L(start) = ceil(5+ (ilength-1)*11);
    start = start+1;
end

for ilength = 1:siz
    L_ind{ilength} = ilength:siz;
end

for icyc = 1:siz
    thisL = L(icyc);
    thisData = dataTrialsCells{icyc};
    x = 1;
    for i = L_ind{icyc}
        thisTrials = dataStructVar.cycNTrials{icyc}
        y = x + thisTrials - 1;
        thisShortData(1:thisL,x:y,:) = thisData(1:thisL,:,:);
        x = y+1;
    end
    cycNoTarget{icyc} = thisShortData;
    clear thisL thisData thisTrials thisShortData x y
end
        
%% responsive cells - baseline vis stim

oneStimResp_data = cycNoTarget{1};
nCells = size(oneStimResp_data,3);
nTrials = size(oneStimResp_data,2);

pmtNoise = max(mean(mean(oneStimResp_data(6:8,:,:),2),3));

for itrial = 1:nTrials
    for icell = 1:nCells
        preStimResp(itrial,icell) = mean(oneStimResp_data(1:5,itrial,icell),1);
    end
end

for itrial = 1:nTrials
    for icell = 1:nCells
        baselineVisStimResp(itrial,icell) = mean(oneStimResp_data(6:end,itrial,icell),1);
    end
end

baselineVisStimResp = bsxfun(@minus, baselineVisStimResp,pmtNoise);

visResp_ttest = ttest(preStimResp,baselineVisStimResp);
visResp_ind = find(visResp_ttest == 1);




%% Sort by and plot all trial lengths after loading experiment


date = '141215';
mouse = 'AW07';
SubNum = '607';
ImgFolder = '004';
time = '1635';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

    % save in
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

edit Load_SBXdataset_fast.m


%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

%register to averaged frames
data_avg = mean(data(:,:,3000:3010),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data, data_avg);
clear data
data = data_reg;
clear data_reg

data = double(data);

    % variables
nTrials = input.trialSinceReset;
nTrials = nTrials-1;
cLeverDown = double(cell2mat(input.cLeverDown));
cLeverDown = cLeverDown(1:end-1);
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cTargetOn = cTargetOn(1:end-1);
cLeverUp = double(cell2mat(input.cLeverUp));
cLeverUp = cLeverUp(1:end-1);
tCyclesOn = double(cell2mat(input.tCyclesOn));
tCyclesOn = tCyclesOn(1:end-1);
% ONms = input.stimOnTimeMs;
% OFFms = input.stimOffTimeMs;
ONfr = input.nFramesOn;
OFFfr = input.nFramesOff;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(1:end-1);
TrialOutcome = input.trialOutcomeCell;
TrialOutcome = TrialOutcome(1:end-1);
Cycles = unique(tCyclesOn);
minCyclesOn = input.minCyclesOn;
maxCyclesOn = input.maxCyclesOn;
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

%avg trial number of frames for trials that were held through target
start = 1;
siz = size(Cycles,2);
for icyc = 1:siz
    cyc_ind = find(tCyclesOn == Cycles(icyc)); 
    thisL = max(absTargetOn(cyc_ind));
    if thisL > 0
        cycTrialL{1,icyc} = thisL+45;
    elseif isnan(thisL)
        thisL = max(absLeverUp(cyc_ind));
        cycTrialL{1,icyc} = thisL+45
    end
    clear thisL cyc_ind
end


    %data cell
for icyc = 1:siz
        start = 1;
        cyc_ind = find(tCyclesOn == Cycles(icyc));
        thisL = cycTrialL{icyc};
        for itrial = cyc_ind
            trial_ind = cLeverDown(itrial)-4:cLeverDown(itrial)+(thisL-5);
            thisCycData(:,:,start:start+thisL-1) = data(:,:,trial_ind);
            cycData{1,icyc} = thisCycData;
            start = start+thisL;
        end
        clear trial_ind
        clear cyc_ind
        clear thisL
        clear thisCycData
end
    %variable cells
 for icyc = 1:siz
     cycNTrials{1,icyc} = size(cycData{1,icyc},3)/cycTrialL{1,icyc};
     
     cyc_ind = find(tCyclesOn == Cycles(icyc));
     thisL = cycTrialL{1,icyc};
     start = 0;
     start2 = 1;
     for itrial = cyc_ind;
        thisLeverDown(1,start2) = absLeverDown + (start*thisL);
        thisTargetOn(1,start2) = absTargetOn(itrial) + (start*thisL);
        thisLeverUp(1,start2) = absLeverUp(itrial) + (start*thisL);
        start = start+1;
        start2 = start2+1;
     end
     cycLeverDown{1,icyc} = thisLeverDown;
     cycTargetOn{1,icyc} = thisTargetOn;
     cycLeverUp{1,icyc} = thisLeverUp;
     cycBlock2ON{1,icyc} = Block2ON(cyc_ind);
     cycTrialOutcome{1,icyc} = TrialOutcome(cyc_ind);
     
     clear cyc_ind thisL thisLeverDown thisTargetOn thisLeverUp 
 end    


%% structure
    for icyc = 1:siz
        dataStruct.cycData(icyc) = cycData(icyc);
        dataStruct.cycTrialL(icyc) = cycTrialL(icyc);
        dataStruct.cycNTrials(icyc) = cycNTrials(icyc);
        dataStruct.cycLeverDown(icyc) = cycLeverDown(icyc);
        dataStruct.cycTargetOn(icyc) = cycTargetOn(icyc);
        dataStruct.cycLeverOn(icyc) = cycLeverUp(icyc);
        dataStruct.cycBlock2ON(icyc) = cycBlock2ON(icyc);
        dataStruct.cycTrialOutcome(icyc) = cycTrialOutcome(icyc);
    end
    
    dataStruct.ONfr = ONfr;
    dataStruct.OFFfr = OFFfr;
%     dataStruct.ONms = ONms;
%     dataStruct.OFFms = OFFms;    
    dataStruct.RateFRperMS = RateFRperMS;
    dataStruct.Cycles = Cycles;
    dataStruct.minCyclesOn = minCyclesOn;
    dataStruct.maxCyclesOn = maxCyclesOn;

date = '141215';
mouse = 'AW07';
ImgFolder = '004';    
dataStruct.mouse = mouse;
dataStruct.date = date;
dataStruct.ImgFolder = ImgFolder;    
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);


%% dF/F

siz = size(dataStruct.Cycles,2);
for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        thisCycData(:,:,itrial) = mean(dataStruct.cycData{icyc}(:,:,(dataStruct.cycTrialL{icyc}.*(itrial-1)+1):(dataStruct.cycTrialL{icyc}.*(itrial-1)+1)+5),3);
    end
    cycF{1,icyc} = thisCycData;
    clear thisCycData
end

for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        start1 = dataStruct.cycTrialL{icyc}.*(itrial-1)+1;
        start2 = dataStruct.cycTrialL{icyc}.*(itrial);
        thisData = dataStruct.cycData{icyc};
        thisDataF = cycF{icyc};
        thisCycData(:,:,start1:start2) = bsxfun(@minus,thisData(:,:,start1:start2),thisDataF(:,:,itrial));
    end
    cycDF{1,icyc} = thisCycData;
    clear thisCycData thisData thisDataF
end

for icyc = 1:siz
    for itrial = 1:dataStruct.cycNTrials{icyc}
        start1 = dataStruct.cycTrialL{icyc}.*(itrial-1)+1;
        start2 = dataStruct.cycTrialL{icyc}.*(itrial);
        thisDataDF = cycDF{icyc};
        thisDataF = cycF{icyc};
        thisCycData(:,:,start1:start2) = bsxfun(@rdivide,thisDataDF(:,:,start1:start2),thisDataF(:,:,itrial));
    end
    cycDFoverF{1,icyc} = thisCycData;
    clear thisCycData thisDataDF thisDataF
end

dataStructDFoverF.cycData = cycDFoverF;
dataStructDFoverF.cycTrialL = dataStruct.cycTrialL; 
dataStructDFoverF.cycNTrials = dataStruct.cycNTrials;
dataStructDFoverF.cycLeverDown = dataStruct.cycLeverDown;
dataStructDFoverF.cycTargetOn = dataStruct.cycTargetOn;
dataStructDFoverF.cycLeverOn = dataStruct.cycLeverOn;
dataStructDFoverF.cycBlock2ON = dataStruct.cycBlock2ON;
dataStructDFoverF.cycTrialOutcome = dataStruct.cycTrialOutcome;
% dataStructDFoverF.ONms = dataStruct.ONms;
% dataStructDFoverF.OFFms = dataStruct.OFFms;
dataStructDFoverF.ONfr = dataStruct.ONfr;
dataStructDFoverF.OFFfr = dataStruct.OFFfr;
dataStructDFoverF.RateFRperMS = dataStruct.RateFRperMS;
dataStructDFoverF.Cycles = dataStruct.Cycles;
dataStructDFoverF.minCyclesOn = dataStruct.minCyclesOn;
dataStructDFoverF.maxCyclesOn = dataStruct.maxCyclesOn;
dataStructDFoverF.mouse = dataStruct.mouse;
dataStructDFoverF.date = dataStruct.date;
dataStructDFoverF.ImgFolder = dataStruct.ImgFolder;


CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataStructDFoverF.mat','dataStructDFoverF');


    
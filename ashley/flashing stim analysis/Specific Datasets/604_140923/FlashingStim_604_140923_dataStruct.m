%% Sort by and plot all trial lengths for i516 141003 Experiment

%cat all datasets
date = '140923';
mouse = 'AW04';
SubNum = '004';
ImgFolder = '005';
time = '1509';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall.tif');
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
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
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

 %% 006
 clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn

date = '140923';
mouse = 'AW04';
SubNum = '004';
ImgFolder = '006';
time = '1515';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall.tif');
data = double(data);

    % variables
nTrials = input.trialSinceReset;
nTrials = nTrials -1;
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
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(1:end-1);
TrialOutcome = input.trialOutcomeCell;
TrialOutcome = TrialOutcome(1:end-1);
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

    %data cell
for icyc = 1:siz
        start = 1;
        cyc_ind = find(tCyclesOn == Cycles(icyc));
        thisL = cycTrialL{icyc};
        for itrial = cyc_ind
            trial_ind = cLeverDown(itrial)-4:cLeverDown(itrial)+(thisL-5);
            thisCycData(:,:,start:start+thisL-1) = data(:,:,trial_ind);
            cycData006{1,icyc} = thisCycData;
            start = start+thisL;
        end
        clear trial_ind
        clear cyc_ind
        clear thisL
        clear thisCycData
end

    %variable cells
 for icyc = 1:siz
     cycNTrials006{1,icyc} = size(cycData006{1,icyc},3)/cycTrialL{1,icyc};
     
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
     cycLeverDown006{1,icyc} = thisLeverDown;
     cycTargetOn006{1,icyc} = thisTargetOn;
     cycLeverUp006{1,icyc} = thisLeverUp;
     cycBlock2ON006{1,icyc} = Block2ON(cyc_ind);
     cycTrialOutcome006{1,icyc} = TrialOutcome(cyc_ind);
     
     clear cyc_ind thisL thisLeverDown thisTargetOn thisLeverUp
 end   

% cat 
for icyc = 1:siz
    tFrames(icyc) = size(cycData{icyc},3);
    thisCycData = cycLeverDown006{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverDown006{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycTargetOn006{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycTargetOn006{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycLeverUp006{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverUp006{icyc} = thisCycData;
    clear thisCycData
end
    
for icyc = 1:siz
    thisCycData = cat(3,cycData{icyc},cycData006{icyc});
    cycData{1,icyc} = thisCycData;
    clear thisCycData
    cycNTrials{1,icyc} = cycNTrials006{icyc}+cycNTrials{icyc};
    cycLeverDown{1,icyc} = cat(2,cycLeverDown{icyc},cycLeverDown006{icyc});
    cycTargetOn{1,icyc} = cat(2,cycTargetOn{1,icyc},cycTargetOn006{1,icyc});
    cycLeverUp{1,icyc} = cat(2,cycLeverUp{1,icyc},cycLeverUp006{1,icyc});
    cycBlock2ON{1,icyc} = cat(2,cycBlock2ON{1,icyc},cycBlock2ON006{1,icyc});
    cycTrialOutcome{1,icyc} = cat(2,cycTrialOutcome{1,icyc},cycTrialOutcome006{1,icyc});
end
 


%% 007
clear cycData006 cycNTrials006 cycLeverDown006 cycTargetOn006 cycLeverUp006 cycBlock2ON006 cycTrialOutcome006 
clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn

date = '140923';
mouse = 'AW04';
SubNum = '004';
ImgFolder = '007';
time = '1527';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall.tif');
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
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(1:end-1);
TrialOutcome = input.trialOutcomeCell;
TrialOutcome = TrialOutcome(1:end-1);
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

    %data cell
for icyc = 1:siz
        start = 1;
        cyc_ind = find(tCyclesOn == Cycles(icyc));
        thisL = cycTrialL{icyc};
        for itrial = cyc_ind
            trial_ind = cLeverDown(itrial)-4:cLeverDown(itrial)+(thisL-5);
            thisCycData(:,:,start:start+thisL-1) = data(:,:,trial_ind);
            cycData007{1,icyc} = thisCycData;
            start = start+thisL;
        end
        clear trial_ind
        clear cyc_ind
        clear thisL
        clear thisCycData
end
    %variable cells
 for icyc = 1:siz
     cycNTrials007{1,icyc} = size(cycData007{1,icyc},3)/cycTrialL{1,icyc};
     
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
     cycLeverDown007{1,icyc} = thisLeverDown;
     cycTargetOn007{1,icyc} = thisTargetOn;
     cycLeverUp007{1,icyc} = thisLeverUp;
     cycBlock2ON007{1,icyc} = Block2ON(cyc_ind);
     cycTrialOutcome007{1,icyc} = TrialOutcome(cyc_ind);
     
     clear cyc_ind thisL thisLeverDown thisTargetOn thisLeverUp
 end 

% cat 
for icyc = 1:siz
    tFrames(icyc) = size(cycData{icyc},3);
    thisCycData = cycLeverDown007{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverDown007{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycTargetOn007{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycTargetOn007{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycLeverUp007{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverUp007{icyc} = thisCycData;
    clear thisCycData
end
    
for icyc = 1:siz
    thisCycData = cat(3,cycData{icyc},cycData007{icyc});
    cycData{1,icyc} = thisCycData;
    clear thisCycData
    cycNTrials{1,icyc} = cycNTrials007{icyc}+cycNTrials{icyc};
    cycLeverDown{1,icyc} = cat(2,cycLeverDown{icyc},cycLeverDown007{icyc});
    cycTargetOn{1,icyc} = cat(2,cycTargetOn{1,icyc},cycTargetOn007{1,icyc});
    cycLeverUp{1,icyc} = cat(2,cycLeverUp{1,icyc},cycLeverUp007{1,icyc});
    cycBlock2ON{1,icyc} = cat(2,cycBlock2ON{1,icyc},cycBlock2ON007{1,icyc});
    cycTrialOutcome{1,icyc} = cat(2,cycTrialOutcome{1,icyc},cycTrialOutcome007{1,icyc});
end
 
clear cycData007 cycNTrials007 cycLeverDown007 cycTargetOn007 cycLeverUp007 cycBlock2ON007 cycTrialOutcome007 
clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn
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
    
    dataStruct.ONms = ONms;
    dataStruct.OFFms = OFFms;
    dataStruct.RateFRperMS = RateFRperMS;
    dataStruct.Cycles = Cycles;
    dataStruct.minCyclesOn = minCyclesOn;
    dataStruct.maxCyclesOn = maxCyclesOn;

date = '140923';
mouse = 'AW04';
ImgFolder = '005+006+007';    
dataStruct.mouse = mouse;
dataStruct.date = date;
dataStruct.ImgFolder = ImgFolder;    
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataStruct.mat','dataStruct');


    
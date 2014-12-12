%% Sort by and plot all trial lengths for i516 141003 Experiment

%cat all datasets
date = '141003';
mouse = '516';
ImgFolder = '002';
time = '1156';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' mouse '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall.tif');
data = double(data);

    % variables
nTrials = input.trialSinceReset;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
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
for icyc = Cycles
    cyc_ind = find(tCyclesOn == Cycles(icyc));
    thisL = max(absTargetOn(cyc_ind));
    if thisL > 0
        cycTrialL{1,icyc} = thisL+40;
    elseif isnan(thisL)
        thisL = max(absLeverUp(cyc_ind));
        cycTrialL{1,icyc} = thisL+40
    end    
    clear thisL cyc_ind
end


    %data cell
for icyc = Cycles
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
 for icyc = Cycles
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

 %% 003
 clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn

date = '141003';
mouse = '516';
ImgFolder = '003';
time = '1159';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' mouse '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall_003.tif');
data = double(data);

    % variables
nTrials = input.trialSinceReset;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

    %data cell
for icyc = Cycles
        start = 1;
        cyc_ind = find(tCyclesOn == Cycles(icyc));
        thisL = cycTrialL{icyc};
        for itrial = cyc_ind
            trial_ind = cLeverDown(itrial)-4:cLeverDown(itrial)+(thisL-5);
            thisCycData(:,:,start:start+thisL-1) = data(:,:,trial_ind);
            cycData003{1,icyc} = thisCycData;
            start = start+thisL;
        end
        clear trial_ind
        clear cyc_ind
        clear thisL
        clear thisCycData
end

    %variable cells
 for icyc = Cycles
     cycNTrials003{1,icyc} = size(cycData003{1,icyc},3)/cycTrialL{1,icyc};
     
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
     cycLeverDown003{1,icyc} = thisLeverDown;
     cycTargetOn003{1,icyc} = thisTargetOn;
     cycLeverUp003{1,icyc} = thisLeverUp;
     cycBlock2ON003{1,icyc} = Block2ON(cyc_ind);
     cycTrialOutcome003{1,icyc} = TrialOutcome(cyc_ind);
     
     clear cyc_ind thisL thisLeverDown thisTargetOn thisLeverUp
 end   

% cat 
for icyc = Cycles
    tFrames(icyc) = size(cycData{icyc},3);
    thisCycData = cycLeverDown003{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverDown003{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycTargetOn003{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycTargetOn003{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycLeverUp003{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverUp003{icyc} = thisCycData;
    clear thisCycData
end
    
for icyc = Cycles
    thisCycData = cat(3,cycData{icyc},cycData003{icyc});
    cycData{1,icyc} = thisCycData;
    clear thisCycData
    cycNTrials{1,icyc} = cycNTrials003{icyc}+cycNTrials{icyc};
    cycLeverDown{1,icyc} = cat(2,cycLeverDown{icyc},cycLeverDown003{icyc});
    cycTargetOn{1,icyc} = cat(2,cycTargetOn{1,icyc},cycTargetOn003{1,icyc});
    cycLeverUp{1,icyc} = cat(2,cycLeverUp{1,icyc},cycLeverUp003{1,icyc});
    cycBlock2ON{1,icyc} = cat(2,cycBlock2ON{1,icyc},cycBlock2ON003{1,icyc});
    cycTrialOutcome{1,icyc} = cat(2,cycTrialOutcome{1,icyc},cycTrialOutcome003{1,icyc});
end
 


%% 004
clear cycData003 cycNTrials003 cycLeverDown003 cycTargetOn003 cycLeverUp003 cycBlock2ON003 cycTrialOutcome003 
clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn

date = '141003';
mouse = '516';
ImgFolder = '004';
time = '1203';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' mouse '-' date '-' time]; 
load (mworks);

    % complete tif
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
data = readtiff('FSall_004.tif');
data = double(data);

    % variables
nTrials = input.trialSinceReset;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

    %data cell
for icyc = Cycles
        start = 1;
        cyc_ind = find(tCyclesOn == Cycles(icyc));
        thisL = cycTrialL{icyc};
        for itrial = cyc_ind
            trial_ind = cLeverDown(itrial)-4:cLeverDown(itrial)+(thisL-5);
            thisCycData(:,:,start:start+thisL-1) = data(:,:,trial_ind);
            cycData004{1,icyc} = thisCycData;
            start = start+thisL;
        end
        clear trial_ind
        clear cyc_ind
        clear thisL
        clear thisCycData
end
    %variable cells
 for icyc = Cycles
     cycNTrials004{1,icyc} = size(cycData004{1,icyc},3)/cycTrialL{1,icyc};
     
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
     cycLeverDown004{1,icyc} = thisLeverDown;
     cycTargetOn004{1,icyc} = thisTargetOn;
     cycLeverUp004{1,icyc} = thisLeverUp;
     cycBlock2ON004{1,icyc} = Block2ON(cyc_ind);
     cycTrialOutcome004{1,icyc} = TrialOutcome(cyc_ind);
     
     clear cyc_ind thisL thisLeverDown thisTargetOn thisLeverUp
 end 

% cat 
for icyc = Cycles
    tFrames(icyc) = size(cycData{icyc},3);
    thisCycData = cycLeverDown004{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverDown004{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycTargetOn004{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycTargetOn004{icyc} = thisCycData;
    clear thisCycData
    thisCycData = cycLeverUp004{icyc};
    thisCycData = bsxfun(@plus,thisCycData,tFrames(icyc));
    cycLeverUp004{icyc} = thisCycData;
    clear thisCycData
end
    
for icyc = Cycles
    thisCycData = cat(3,cycData{icyc},cycData004{icyc});
    cycData{1,icyc} = thisCycData;
    clear thisCycData
    cycNTrials{1,icyc} = cycNTrials004{icyc}+cycNTrials{icyc};
    cycLeverDown{1,icyc} = cat(2,cycLeverDown{icyc},cycLeverDown004{icyc});
    cycTargetOn{1,icyc} = cat(2,cycTargetOn{1,icyc},cycTargetOn004{1,icyc});
    cycLeverUp{1,icyc} = cat(2,cycLeverUp{1,icyc},cycLeverUp004{1,icyc});
    cycBlock2ON{1,icyc} = cat(2,cycBlock2ON{1,icyc},cycBlock2ON004{1,icyc});
    cycTrialOutcome{1,icyc} = cat(2,cycTrialOutcome{1,icyc},cycTrialOutcome004{1,icyc});
end
 
clear cycData004 cycNTrials004 cycLeverDown004 cycTargetOn004 cycLeverUp004 cycBlock2ON004 cycTrialOutcome004 
clear data nTrials cLeverDown cTargetOn cLeverUp tCyclesOn Block2ON TrialOutcome absLeverDown absLeverUp absTargetOn
%% structure
    for icyc = Cycles
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

date = '141003';
mouse = '516';
ImgFolder = '002+003+004';    
dataStruct.mouse = mouse;
dataStruct.date = date;
dataStruct.ImgFolder = ImgFolder;    
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataStruct.mat','dataStruct');


    
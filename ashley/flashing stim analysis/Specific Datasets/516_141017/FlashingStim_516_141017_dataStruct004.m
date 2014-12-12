%% Sort by and plot all trial lengths for i516 141003 Experiment

%one 27000 frame dataset
date = '141017';
mouse = '516';
ImgFolder = '004';
time = '1026';

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

date = '141017';
mouse = '516';
ImgFolder = '004';    
dataStruct.mouse = mouse;
dataStruct.date = date;
dataStruct.ImgFolder = ImgFolder;    
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
save('dataStruct.mat','dataStruct');


    
mouse = 'AW04';
SubNum = '604';
date = '140923';
ImgFolder = '005+006+007';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat')

%get average response to vis stim for each cycle

%% dirTuning dataset
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

    %avg all trials all cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    for itrial = 1:dataStructVar.cycNTrials{icyc}
        thisDataAvgCells(:,itrial) = mean(thisData(:,itrial,:),3);
    end
    thisDataAvgTrials = mean(thisDataAvgCells,2);
    cycAvgTC{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells thisDataAvgTrials
end

    %FIG for traces
nCycles = max(dataStructVar.Cycles);
figure;
color = colormap(hsv(nCycles));
for icyc = 1:siz
    cycColor = color(icyc,:);
    plot(cycAvgTC{icyc},'Color',cycColor);
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    hold on
end

%find significant responses to vis stim, paired t-test to test visual response
for icyc = 1:siz
    thisL = dataStructDFoverF.cycTrialL{icyc};
    trials = dataStructDFoverF.cycNTrials{icyc};
    thisTC = dataTrialsCells{icyc};
    for itrial = 1:trials
        for icell = 1:nCells
            this_last(:,itrial,icell) = mean(thisTC(end-19:end,itrial,icell),1);
        end
    end
    dataTrialsCells_last{icyc} = this_last;
    clear thisL trials thisTC this_last
end

for icyc = 1:siz
    thisL = dataStructDFoverF.cycTrialL{icyc};
    trials = dataStructDFoverF.cycNTrials{icyc};
    thisTC = dataTrialsCells{icyc};
    for itrial = 1:trials
        for icell = 1:nCells
            this_first(:,itrial,icell) = mean(thisTC(1:5,itrial,icell),1);
        end
    end
    dataTrialsCells_first{icyc} = this_first;
    clear thisL trials thisTC this_first
end

    %paired t-test
for icyc = 1:siz
    this_first = dataTrialsCells_first{icyc};
    this_last = dataTrialsCells_last{icyc};
    first2last_ttest{icyc} = ttest(this_first,this_last);
    visualresp_ind{icyc} = find(first2last_ttest{icyc} == 1);
end


    %avg all trials, visually responsive cells
for icyc = 1:siz
    thisData = dataTrialsCells{icyc};
    cells  = visualresp_ind{icyc};
    for itrial = 1:dataStructDFoverF.cycNTrials{icyc}
        thisDataAvgCells_Rsp(:,itrial) = mean(thisData(:,itrial,cells),3);
    end
    thisDataAvgTrials = mean(thisDataAvgCells_Rsp,2);
    cycAvgTC_Rsp{1,icyc} = thisDataAvgTrials;
    clear thisData thisDataAvgCells_Rsp thisDataAvgTrials_Rsp cells
end

    %FIG for traces
nCycles = max(dataStructDFoverF.Cycles);
figure;
color = colormap(hsv(nCycles));
for icyc = 1:siz
    cycColor = color(icyc,:);
    plot(cycAvgTC_Rsp{icyc},'Color',cycColor);
    hold on
    for i = 1:nCycles+1
    v = 5+ (i-1)*10.5;
    vline(v,'k');
    end
    hold on
end


mouse = '516';
date = '141003';
ImgFolder = '002+003+004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat')
load('dataMasks.mat')

%% Save behavior and imaging parameters sorted by cycles
dataStructVar.cycTrialL = dataStructDFoverF.cycTrialL;
dataStructVar.cycNTrials = dataStructDFoverF.cycNTrials;
dataStructVar.cycLeverDown = dataStructDFoverF.cycLeverDown;
dataStructVar.cycTargetOn = dataStructDFoverF.cycTargetOn;
dataStructVar.cycLeverUp = dataStructDFoverF.cycLeverOn;
dataStructVar.cycBlock2ON = dataStructDFoverF.cycBlock2ON;
dataStructVar.cycTrialOutcome = dataStructDFoverF.cycTrialOutcome;
dataStructVar.ONms = dataStructDFoverF.ONms;
dataStructVar.OFFms = dataStructDFoverF.OFFms;
dataStructVar.RateFRperMS = dataStructDFoverF.RateFRperMS;
dataStructVar.Cycles = dataStructDFoverF.Cycles;
dataStructVar.minCyclesOn = dataStructDFoverF.minCyclesOn;
dataStructVar.maxCyclesOn = dataStructDFoverF.maxCyclesOn;
dataStructVar.mouse = dataStructDFoverF.mouse;
dataStructVar.date = dataStructDFoverF.date;
dataStructVar.ImgFolder = dataStructDFoverF.ImgFolder;
dataStructVar.maxDFoverF = dataStructDFoverF.maxDFoverF;

save('dataStructVar.mat','dataStructVar');

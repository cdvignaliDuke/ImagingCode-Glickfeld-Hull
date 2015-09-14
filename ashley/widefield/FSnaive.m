mouse = 'AW36';
subN = '636';
date = '150911';
time = '1420';
experimentname = 'msFRONT_100exp_10Hz_FSAVnaive_s3_1';

cd('Y:\home\andrew\Behavior\Data')
load(['data-i' subN '-' date '-' time]);

cd(['Z:\data\AW36\widefield imaging\' date '_' mouse '\' experimentname]);
data = double(readtiff([experimentname '_MMStack.ome.tif']));


%% find dF/F for entire frame - only one cycle used
trialLength = 20+(input.maxCyclesOn*(input.nFramesOff+input.nFramesOn));
cLeverDown = cell2mat_padded(input.cLeverDown);
imgrate = input.frameRateHz;


dataDF = zeros(size(data,1),size(data,2),input.trialSinceReset,trialLength);
dataDFoverF = zeros(size(data,1),size(data,2),input.trialSinceReset,trialLength);

for itrial = 1:input.trialSinceReset
    
    dataDF(:,:,itrial,:) = bsxfun(@minus, data(:,:,cLeverDown(itrial)-9:cLeverDown(itrial)+trialLength-10),mean(data(:,:,cLeverDown(itrial)-9:cLeverDown(itrial)),3));
    dataDFoverF(:,:,itrial,:) = bsxfun(@rdivide, dataDF(:,:,itrial,:), mean(data(:,:,cLeverDown(itrial)-9:cLeverDown(itrial)),3));
end

%%
% avgDFoverFimg = squeeze(mean(dataDFoverF,3));
% writetiff(avgDFoverFimg,'avgAllTrialsFullFrame');

%% choose ROIs
max_dF = max(max(dataDF,[],4),[],3);
max_dFoverF = max(max(dataDFoverF,[],4),[],3);
figure; imagesq(max_dF); colormap(gray)

% save max dF/F image
% writetiff(max_dFoverF,'maxDFoverFAllTrialsFullFrame');

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);

dataTC = stackGetTimeCourses(data,mask_cell);
figure; tcOffsetPlot(dataTC)

%% dF/F for timecourses

dataTC_DF = zeros(trialLength,size(dataTC,2),input.trialSinceReset);
dataTC_DFoverF = zeros(trialLength,size(dataTC,2),input.trialSinceReset);

for itrial = 1:input.trialSinceReset
    dataTC_DF(:,:,itrial) = bsxfun(@minus, dataTC(cLeverDown(itrial)-9:(cLeverDown(itrial)+trialLength-10),:),mean(dataTC(cLeverDown(itrial)-9:cLeverDown(itrial),:),1));
    dataTC_DFoverF(:,:,itrial) = bsxfun(@rdivide, dataTC_DF(:,:,itrial),mean(dataTC(cLeverDown(itrial)-9:cLeverDown(itrial),:),1));
end

%% plot timecourses
figure;
plot(mean(dataTC_DFoverF,3))
vline(10, 'k')
for i = 1:input.maxCyclesOn-1
    vline(10+((input.nFramesOff+input.nFramesOn)*i),'k:')
end
vline(10+((input.nFramesOff+input.nFramesOn)*input.maxCyclesOn),'c')
title([mouse ' Naive FS-AV - ' num2str(imgrate) 'Hz, avg all trials']);
xlabel('frames')
ylabel('dF/F')

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

filesave = 'Z:\Analysis\AW36\widefield imaging\150911';

print([filesave ['\' [mouse ' avg all trials, rough ROIs,10Hz_new' date] '.pdf']], '-dpdf')    

%% plot aud vs vis trials timecourses

block2 = cell2mat_padded(input.tBlock2TrialNumber);
visInd = find(block2 == 0);
audInd = find(block2 == 1);

figure;
plot(mean(mean(dataTC_DFoverF(:,:,visInd),3),2),'g')
hold on
plot(mean(mean(dataTC_DFoverF(:,:,audInd),3),2),'k')
vline(10, 'k')
for i = 1:input.maxCyclesOn-1
    vline(10+((input.nFramesOff+input.nFramesOn)*i),'k:')
end
vline(10+((input.nFramesOff+input.nFramesOn)*input.maxCyclesOn),'c')
title([mouse ' Naive FS-AV - ' num2str(imgrate) 'Hz, v = gr, a = bl trials']);
xlabel('frames')
ylabel('dF/F')

filesave = 'Z:\Analysis\AW36\widefield imaging\150911';

print([filesave ['\' [mouse ' avg all trials, avg ROIs,10Hz_new,VvsA' date] '.pdf']], '-dpdf')    



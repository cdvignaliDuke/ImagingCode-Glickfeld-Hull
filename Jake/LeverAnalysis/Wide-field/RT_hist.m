%RT_Histogram    Plots histograms of reaction times relative to cue onset
%for all animals in days. segregates accoring to long vs short randreq hold
%times. 
clear;
%all 16 animals. No bx criteria
days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41'}; %'150718_img27', '150719_img27',
%days = {'160515_img48', '160516_img47', '150718_img27'}; %naive mice    , '150716_img27'
BEHAVE_DIR = 'Z:\Data\WidefieldImaging\GCaMP\behavior\';
BX_OUTPUT_DIR = 'Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\';
reactTimesCellShort = struct(); 
reactTimesCellLong = struct(); 
trialNumsShort = [];
trialNumsLong = [];
daysShort = [];
daysLong = [];
%open each bx day file, determine if long or short day. Plot hist for that
%day. Store variables to be used for the averages. 
for kk=1:length(days)
   day_name  =  days{kk};
   bfile = dir([BEHAVE_DIR 'data-*i9' days{kk}(end-1:end) '-' days{kk}(1:6) '*' ]);
   behave_dest = [BEHAVE_DIR bfile.name];
   assert(length(bfile)) =1;
   b_data = load(behave_dest);
   %disp([days{kk}, ' ', num2str(b_data.input.randReqHoldMaxMs), ' ', num2str(b_data.input.fixedReqHoldTimeMs)])
   load([BX_OUTPUT_DIR days{kk} '_bx_outputs.mat'], 'frame_info')
   f_trial = frame_info.f_frame_trial_num+1;
   l_trial = frame_info.l_frame_trial_num-1;
   reactTimesMs = double(cell2mat(b_data.input.reactTimesMs(f_trial:l_trial)));
   holdTimesMs = double(cell2mat(b_data.input.holdTimesMs(f_trial:l_trial)));
   if b_data.input.randReqHoldMaxMs < 2000;
       xbins = [-1500:50:1000];
       xbinsHold = [0:50:2500];
       daysShort = [daysShort, days(kk)];
   else
       xbins = [-5000:50:1000];
       xbinsHold = [0:50:6000];
       daysLong = [daysLong, days(kk)];
   end
   countsHold = hist(holdTimesMs, xbinsHold);
   countsHold = countsHold/sum(countsHold);
   counts = hist(reactTimesMs, xbins);
   counts = counts/sum(counts);
%    figure; subplot(1,2,1); plot(xbinsHold,cumsum(countsHold))
%    title([days{kk}(1:6), ' ', days{kk}(end-4:end), ' 50ms bins n=', num2str(length(reactTimesMs))]);
%    xlabel('Hold Time (Ms)');
%    ylabel('fraction of trials');
%    ylim([0 1])
%    if b_data.input.randReqHoldMaxMs < 2000;
%        xlim([0 2500])
%        else
%        xlim([0 6000]) 
%    end
   %figure; hist(reactTimesMs,xbins); 
   figure;
%    subplot(1,2,2); plot(xbins,cumsum(counts))
%    title([days{kk}(1:6), ' ', days{kk}(end-4:end), ' 50ms bins n=', num2str(length(reactTimesMs))]);
%    xlabel('Reaction Time (Ms)');
%    ylabel('fraction of trials');
%    ylim([0 1])
%    vline(0);
%    vline(500);
%    if b_data.input.randReqHoldMaxMs < 2000;
%        xlim([-1500 1000])
%        reactTimesCellShort.(['i', daysShort{end}])=reactTimesMs;    %collect react times for population histogram   
%        trialNumsShort = [trialNumsShort, length(reactTimesMs)];    %collect trial number for population histogram
%    else
%        xlim([-5000 1000])
%        reactTimesCellLong.(['i', daysLong{end}])=reactTimesMs;      
%        trialNumsLong = [trialNumsLong, length(reactTimesMs)]; 
%    end
   plot([xbins(1:end-1)+25], diff(counts));
   title([days{kk}(1:6), ' ', days{kk}(end-4:end), ' 50ms bins n=', num2str(length(reactTimesMs))]);
   xlabel('Reaction Time (Ms)');
   ylabel('slope of line');
   ylim([0 0.4])
   vline(0);
   vline(500);
end
clear day_name bfile behave_dest b_data f_trial l_trial reactTimesMs frame_info assert; 

%short randHold days
minTrialNumShort = min(trialNumsShort);
clippedReactTimes = NaN(length(daysShort),minTrialNumShort); 
for kk=1:length(daysShort);
    if length(reactTimesCellShort.(['i',daysShort{kk}])) > minTrialNumShort;
        clippedReactTimes(kk,:) = randsample(reactTimesCellShort.(['i',daysShort{kk}]),minTrialNumShort);
    else
        clippedReactTimes(kk,:) = reactTimesCellShort.(['i',daysShort{kk}]);
    end
end
reshapedReactTimes = reshape(clippedReactTimes, 1, minTrialNumShort*length(daysShort));
%countsShort = hist(reshapedReactionTimes, [-1500:50:1000]);
figure; hist(reshapedReactTimes,[-1500:50:1000]);
title(['Population Reaction Times randHold=1000 ', num2str(minTrialNumShort), 'trials per animal']);
xlabel('Reaction Time (Ms)');
ylabel('number of trials per bin');
xlim([-1500 1000]);
vline(0);

%long randHold days
minTrialNumLong = min(trialNumsLong);
clippedReactTimes = NaN(length(daysLong),minTrialNumLong); 
for kk=1:length(daysLong);
    if length(reactTimesCellLong.(['i',daysLong{kk}])) > minTrialNumLong;
        clippedReactTimes(kk,:) = randsample(reactTimesCellLong.(['i',daysLong{kk}]),minTrialNumLong);
    else
        clippedReactTimes(kk,:) = reactTimesCellLong.(['i',daysLong{kk}]);
    end
end
reshapedReactTimes = reshape(clippedReactTimes, 1, minTrialNumLong*length(daysLong));
figure; hist(reshapedReactTimes,[-5000:50:1000]);
title(['Population Reaction Times randHold=4500 ', num2str(minTrialNumLong), 'trials per animal']);
xlabel('Reaction Time (Ms)');
ylabel('number of trials per bin');
xlim([-5000 1000]);
vline(0);
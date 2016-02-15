bxSumFig = figure;
%% concatenate input structures
mworksLocation = 'Y:\home\andrew\Behavior\Data';

for ifile = 1:size(taskTimeAll,1)
            if ifile == 1
                input_temp = mwLoadData(fullfile(mworksLocation, ['data-' 'i' SubNum '-' date '-' taskTimeAll(ifile,:) '.mat']), [], []);
            else
                next_input = mwLoadData(fullfile(mworksLocation, ['data-' 'i' SubNum '-' date '-' taskTimeAll(ifile,:) '.mat']), [], []);
                fn1 = fieldnames(input_temp);
                fn2 = fieldnames(next_input);
                missingVariables = [];
                if length(fn2) < length(fn1) %% in case input structures have different number of field names
                    missingVariables = setdiff(fn1,fn2);
                    for ivar = 1:length(missingVariables)
                        next_input.(char(missingVariables{i})) = [];
                    end 
                    input_temp = [input_temp next_input];
                elseif length(fn2) > length(fn1)
                    missingVariables = setdiff(fn2 , fn1);
                    for ivar = 1:length(missingVariables)
                        input_temp.(char(missingVariables{ivar})) = [];
                    end
                    input_temp = [input_temp next_input];
                else
                    input_temp = [input_temp next_input];
                end
            end
            
end
input_temp = concatenateDataBlocks(input_temp);

%% trial outcomes and target values
failureIx = strcmp(input_temp.trialOutcomeCell, 'failure');
missedIx = strcmp(input_temp.trialOutcomeCell, 'ignore');
successIx = strcmp(input_temp.trialOutcomeCell, 'success');
pctEarly = sum(failureIx,2)./length(failureIx);
avgHR = sum(successIx)/(sum(successIx) + sum(missedIx));
pctMiss = sum(missedIx)/length(missedIx);

gratingDirectionDeg = celleqel2mat_padded(input_temp.tGratingDirectionDeg);

oris = unique(gratingDirectionDeg);



%% calc hit rates
HR_V = zeros(1,length(oris)-1);
HR_AV = [];
for i = 1:length(oris)
   ind = find(gratingDirectionDeg == oris(i));
   if i == 1
       HR_AV = sum(successIx(ind)/(sum(successIx(ind))+sum(missedIx(ind))));
   else
       HR_V(i-1) = sum(successIx(ind)/(sum(successIx(ind))+sum(missedIx(ind))));
   end
end


%% plot hit rates
figure(bxSumFig);
subplot(2,1,1);
xAxis_V = oris(2:end);
xAxis_AV = 100;
% try 
%     xAxis_AV = input_temp.block2SoundTargetAmplitude;
% catch
%     xAxis_AV = input_temp.soundTargetAmplitude*100;
% end

ax = gca;
plot(xAxis_V,HR_V,'g.','MarkerSize',30)
hold on
plot(xAxis_AV,HR_AV,'k.','MarkerSize',30)
hold on
set(gca, 'xscale', 'log');
xlim([0 100])
ylim([0 1.1])
xlabel('Orientation change (deg)')
ylabel('Hit rate')
xtick = oris;
set(ax,'XTick',xtick)
title({[mouse '-' date];'hit rate summary'})
legend({'vis';'aud'},'Location','southeastoutside')

%% calc cumulative FA, miss, and correct
lag = 10;
cumFA = tsmovavg(cumsum(failureIx)./(1:length(failureIx)),'s',lag,2);
cumMiss = tsmovavg(cumsum(missedIx)./(1:length(failureIx)),'s',lag,2);
cumSuccess = tsmovavg(cumsum(successIx)./(1:length(failureIx)),'s',lag,2);

%% calc cumulative FA, miss, and correct
figure(bxSumFig)
subplot(2,1,2);
ax = gca;
plot(cumFA,'c')
hold on
plot(cumMiss,'m')
hold on
plot(cumSuccess,'k')
ylim([0 1])
xtick = cumsum(input_temp.trialsSinceReset);
xlabel('trial number')
set(ax,'XTick',xtick)
ylabel('smoothed percent')
title({[mouse '-' date];'cumulative summary'})
legend({'FA';'miss';'hit'}, 'Location', 'southeastoutside')

%% save fig
figure(bxSumFig)
print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'bxSummary'), '-dpdf');

%% if there are catch trials...

%% if there are multiple sound target amplitudes
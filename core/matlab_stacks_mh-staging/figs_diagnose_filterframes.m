function figs_diagnose_filterframes(stack, bwMask, frameTimeMs, interpFrames, ...
                                    stimTs)
%
%  draw some figures to make sure you're using the right interp frames
%
%$Id: figs_diagnose_filterframes.m 78 2008-01-10 15:38:51Z histed $

%% process options
ts = stimTs;

%% process tcourses
tc = stack_get_tcourses(stack, bwMask, 'mean');
[nCells, nFrames] = size(tc);
filtOpts = { 'FrameTimeMs', frameTimeMs, ...
             'BaselineFrames', interpFrames, ...
             'DoHighpassFilter', true, ...
           };

% no filt
tcFN = stack_filter_tcourses(tc, filtOpts{1:4}, ...
                             'DoHighpassFilter', false);
tcMat = vec2padded(tcFN, ts.blockStartNs, [], ts.blockEndNs, NaN);



% filt butter - full
tcFB = stack_filter_tcourses(tc, filtOpts{:}, ...
                             'FilterType', 'butter');
tcMatB = vec2padded(tcFB, ts.blockStartNs, [], ts.blockEndNs, NaN);

% filt interp
tcFI = stack_filter_tcourses(tc, filtOpts{:}, ...
                             'FilterType', 'interp', ...
                             'FilterInterpFrames', interpFrames);
tcMatI = vec2padded(tcFI, ts.blockStartNs, [], ts.blockEndNs, NaN);


figH = figure;

spH = subplot(2,2,1);
hold on;
plot(tcMatB);
plot(mean(tcMatB',1), ...
     'LineWidth', 3, ...
     'Color', 'k');
set(gca, 'XLim', [1 ts.nFramesInBlock]);
title('Butterworth');

spH = subplot(2,2,2);
hold on;
plot(tcMatI);
plot(mean(tcMatI',1), ...
     'LineWidth', 3, ...
     'Color', 'k');
set(gca, 'XLim', [1 ts.nFramesInBlock]);
title('interp filtering');

spH = subplot(2,2,3);
hold on;
plot(tcMat);
plot(mean(tcMat',1), ...
     'LineWidth', 3, ...
     'Color', 'k');
set(gca, 'XLim', [1 ts.nFramesInBlock]);
title('no filtering');

spH = subplot(2,2,4);
hold on;
noM = mean(tcMat');
interM = mean(tcMatI');
butterM = mean(tcMatB');
pH(1) = plot(noM, 'rx');
pH(2) = plot(interM, 'bx');
pH(3) = plot(butterM, 'gx');
pH(4) = plot(interM-noM, 'kx-');

legH = legend(pH, {'none', 'interp', 'full', 'interp-none'});

set(gca, 'XLim', [1 ts.nFramesInBlock]);
title('means and (interp - no filter) diff');

suptitle2(make_title(': any mean shifts induced by filtering?'));


%% new scatter figure
mI = mean(tcMatI',1);
mB = mean(tcMatB',1);
sI = std(tcMatI');
sB = std(tcMatB');
figH = figure;
hold on;
iPH = plot(mI,sI, 'bo');
bPH = plot(mB,sB, 'ro');
xlabel('mean of frame');
ylabel('std of frame');
legH = legend([iPH bPH], {'Interp', 'Butter(full)'});

% eb
errorbar(mean(mI), mean(sI), std(sI), 'b');
ebH = errorbar(mean(mB), mean(sB), std(sB), 'r');

title(make_title(': have you used enough interp frames to avoid adding noise?'));

fprintf(1, '%s: std mean over frames: interp %6.4g, full(butter) %6.4g\n', ...
        mfilename, mean(sI), mean(sB));


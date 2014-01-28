%% get a sample stack
cd ../histed_estim3
[stack cs] = cat061012_get_stack('estim18');
cd ../estim_iterator


%%%%%%%%%%%%%%%%
%% various consts
uo.StimEvery = cs.stimEvery;

[nRows nCols nFrames] = size(stack);
% stim averaging consts
preStimNs = floor(uo.StimEvery/2);
postStimNs = ceil(uo.StimEvery/2);
tsAvg = stim_make_frame_ns(nFrames, uo.StimEvery, ...
                           preStimNs, postStimNs, 1);


%% test movie
[movFr scaleRange] = make_stim_avgmovie(stack, ...
                                        'FrameNsToAverage', tsAvg.frameNsToAverage, ...
                                        'BlockStimNs', tsAvg.blockStimNs, ...
                                        'DoPlot', true, ...
                                        'MovieType', 'dFOF', ...
                                        'DFOFBlockBaseNs', 1:5, ...
                                        'PlotColormap', hot(256), ...
                                        'PlotDoColorbar', true);

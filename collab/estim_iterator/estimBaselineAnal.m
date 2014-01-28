function estimBaselineAnal(stack, varargin)

defs = { 'OutputDir', [], ...
         'ExptName', [], ...
         'SeriesName', [], ...
         'DoTruncate', [], ...
         'StimEvery', [], ...
         'Current', [], ...
         'FrameTimeMs', [], ...
         'TrainTimeMs', [], ...
         'NStimsInSer', [], ...
       };
         
uo = stropt2struct(stropt_defaults(defs, varargin));


% brilliant.  Now iterate over several things to do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
%% stack size
[nRows nCols nFrames] = size(stack);

%%%%%%%%%%%%%%%%
%% deal with truncated stacks
if isnan(uo.DoTruncate) ...   % if missing, try it anyway
        || uo.DoTruncate == true  % and if specified, definitely do it
    maxStims = floor((nFrames-(uo.StimEvery-1)) ./ uo.StimEvery);
    maxReps = floor(maxStims ./ uo.NStimsInSer);
    assert(maxReps >= 1, ...
           'Too few frames for even 1 rep');
    nStimsToKeep = maxReps*uo.NStimsInSer;
    nFramesToKeep = uo.StimEvery*nStimsToKeep+(uo.StimEvery-1);
    assert(nFramesToKeep <= nFrames, 'bug');
    if nFramesToKeep < nFrames;
        assert(nFramesToKeep > nFrames/2, ...
               'Truncated by a factor of 2?');
        stack = stack(:,:,1:nFramesToKeep);
        nFrames = nFramesToKeep;
        fprintf(1, '\n    *** Truncated to %d frames from %d\n', ...
                nFramesToKeep, nFrames);
    end
else
    fprintf(1, '\n    ***   Skipping truncate\n');
end

    




%%%%%%%%%%%%%%%%
%% stim averaging consts
preStimNs = floor(uo.StimEvery/2);
postStimNs = ceil(uo.StimEvery/2);
tsAvg = stim_make_frame_ns(nFrames, uo.StimEvery, ...
                           preStimNs, postStimNs, uo.NStimsInSer);
tsAll = stim_make_frame_ns(nFrames, uo.StimEvery, ...
                           preStimNs, postStimNs, 1);
tsInterp = stim_make_frame_ns(nFrames, uo.StimEvery, ...
                              preStimNs, 0, 1);
interpFrNs = sort(reshape(tsInterp.frameNsToAverage(:,1:end-1), 1, []));

%%%%%%%%%%%%%%%%
%% fov
fov = mean(stack,3);
fovS = imScale(fov, [], [0 255], 'uint8');
fovL = stack_localcontrastadj(fov, 31, 5);
fovLS = imScale(fovL, [], [0 255], 'uint8');
imwrite(fovS, fullfile(uo.OutputDir, 'fov.png'));
imwrite(fovLS, fullfile(uo.OutputDir, 'fovLocal.png'));

%%%%%%%%%%%%%%%%
%% overall pix timecourse
pM = squeeze(mean(mean(stack,2),1));
xvals = (1:nFrames)*(uo.FrameTimeMs/1000);
figH = figure;
hold on;
plot(xvals, pM, '-x');
vert_lines(xvals(tsAll.stimNs));
title(sprintf('%s: %s %s, frameTime %dms', ...
              'all pix mean', uo.ExptName, uo.SeriesName, uo.FrameTimeMs));
xlabel('Time (s)');
ylabel('raw F');

set(figH, 'Visible', 'off');
outName = fullfile(uo.OutputDir, 'all_pix_mean');
exportfigPrint(figH, outName, ...
               'FileFormat', 'pdf', ...
               'Size', 10*[1 0.75]);
close(figH);


%%%%%%%%%%%%%%%%
%% overall pix AVERAGE timecourse
pMF = tcFilter(pM, ...
               'FrameTimeMs', uo.FrameTimeMs, ...
               'DoHighpassFilter', true, ...
               'FilterType', 'interp', ...
               'BaselineFrames', interpFrNs, ...
               'FilterInterpFrames', interpFrNs, ...
               'DoDFOF', true, ...
               'FrameNsToAverage', tsAvg.frameNsToAverage, ...
               'DoPostSmooth', false);


xvals = ( (1:tsAvg.nFramesInBlock) - tsAvg.blockStimNs(1)).*(uo.FrameTimeMs/1000);
figH = figure;
hold on;
plot(xvals, pMF, '-x');
vert_lines(xvals(tsAvg.blockStimNs));
title(sprintf('%s: %s %s', ...
              'all pix average mean', uo.ExptName, uo.SeriesName));
xlabel('Time (s)');
ylabel('dF/F_0');

set(figH, 'Visible', 'off');
outName = fullfile(uo.OutputDir, 'all_pix_avg_mean');
exportfigPrint(figH, outName, ...
               'FileFormat', 'pdf', ...
               'Size', 10*[1 0.75]);
close(figH);






%%%%%%%%%%%%%%%%
%% stim_pixmap
for iS=1:uo.NStimsInSer
    tCurr = uo.Current(iS);
    assert(length(uo.Current) == uo.NStimsInSer, 'listfile bug');
    % one for each stim current
    [crap smOutS figH] ...
        = make_stim_pixmap(stack, ...
                           'StimFrNs', tsAvg.stimNs(iS), ...
                           'NBaseFrames', 5, ...
                           'DoPlot', true, ...
                           'CLim', [-0.1 0.35], ...
                           'DoMapSmooth', true, ...
                           'MapSmoothSigma', 3, ...
                           'DoDFOF', true, ...
                           'PrintStatus', false);
    title(sprintf('\\bf%s\\rm: %s %s, %duA (curr %d/%d), train %dms', ...
                  'make\_stim\_pixmap', uo.ExptName, uo.SeriesName, ...
                  tCurr, iS, uo.NStimsInSer, ...
                  uo.TrainTimeMs));

    set(figH, 'Visible', 'off');
    outName = fullfile(uo.OutputDir, ...
                       sprintf('make_stim_pixmap-curr%02d', iS));
    exportfigPrint(figH, outName, ...
                   'FileFormat', 'pdf', ...
                   'Size', 10*[1 0.75]);
    close(figH);
end

%%%%%%%%%%%%%%%%
%% PST movie

% dF
[movFr scaleRangeDF] ...
    = make_stim_avgmovie(stack, ...
                         'FrameNsToAverage', tsAvg.frameNsToAverage, ...
                         'BlockStimNs', tsAvg.blockStimNs, ...
                         'DoPlot', false, ...
                         'MovieType', 'dF');
outName = fullfile(uo.OutputDir, 'stim_avgmovie-dF');
stackWriteAvi(movFr, outName, 5, gray(256));

% dFOF
[movFr scaleRangeDFOF] ...
    = make_stim_avgmovie(stack, ...
                         'FrameNsToAverage', tsAvg.frameNsToAverage, ...
                         'BlockStimNs', tsAvg.blockStimNs, ...
                         'DoPlot', false, ...
                         'MovieType', 'dFOF', ...
                         'DFOFBlockBaseNs', 1:5);
outName = fullfile(uo.OutputDir, 'stim_avgmovie-dFOF');
stackWriteAvi(movFr, outName, 5, hot(256));

%%%%%%%%%%%%%%%%
%% runinfo.txt
% save important data
outTxtName = fullfile(uo.OutputDir, 'runinfo.txt');
fid = fopen(outTxtName, 'w+');
fprintf(fid, '\n*** %s %s\n', mfilename, datestr(now));
fprintf(fid, 'Expt: %s %s\n', uo.ExptName, uo.SeriesName);
fprintf(fid, 'nRows, nCols, nFrames: %s\nstimEvery %3d\n', ...
        mat2str(size(stack)), uo.StimEvery); 
fprintf(fid, 'Listfile: Current %suA, TrainTimeMs %d, FrameTimeMs %d\n', ...
        mat2str(uo.Current(:)'), uo.TrainTimeMs, uo.FrameTimeMs);
fprintf(fid, 'Movie scale ranges: dF %s;   dF/F: %s\n', ...
        mat2str(scaleRangeDF(:)'), mat2str(scaleRangeDFOF(:)'));
fclose(fid);


                                 





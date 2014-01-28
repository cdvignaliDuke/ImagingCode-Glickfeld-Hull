function estimBaselineAnal2(listfileEntry, outputDir, doExport, mapSmoothSigma, pixmapCLim)

if nargin < 3, doExport = true; end
if nargin < 4, mapSmoothSigma = 3; end
if nargin < 5, pixmapCLim = [-0.1 0.35]; end

le = listfileEntry;         

% constants
%mapSmoothSigma = 5; %10;
doMovie = false;  

% brilliant.  Now iterate over several things to do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filePrefix = strcat(le.ExptName, '_', le.SeriesName, '-');

%%%%%%%%%%%%%%%%
%% le.greenStack size
[nRows nCols nFrames] = size(le.greenStack);

%%%%%%%%%%%%%%%%
%% stim averaging consts
as = estimSubAvgConsts(le.StimEvery, nFrames, le.NStimsInSer);

%%%%%%%%%%%%%%%%
%% fov
fov = mean(le.greenStack,3);
fovS = imScale(fov, [], [0 255], 'uint8');
fovL = stack_localcontrastadj(fov, 31, 5);
fovLS = imScale(fovL, [], [0 255], 'uint8');
if doExport
    imwrite(fovS, fullfile(outputDir, [filePrefix 'fov.png']));
    imwrite(fovLS, fullfile(outputDir, [filePrefix 'fovLocal.png']));
else
    % red for display
    figH = figure; imagesc(le.fovRGB); axis square;
    figH = figure; imagesc(le.fovN_RGB); axis square;    
end


%%%%%%%%%%%%%%%%
%% overall pix timecourse
pM = squeeze(mean(mean(le.greenStack,2),1));
xvals = (1:nFrames)*(le.FrameTimeMs/1000);
figH = figure;
hold on;
plot(xvals, pM, '-x');
vert_lines(xvals(as.tsOne.stimNs));
title(sprintf('%s: %s %s, frameTime %dms', ...
              'all pix mean', le.ExptName, le.SeriesName, le.FrameTimeMs));
xlabel('Time (s)');
ylabel('raw F');

if doExport
    set(figH, 'Visible', 'off');
    outName = fullfile(outputDir, [filePrefix 'all_pix_mean']);
    exportfigPrint(figH, outName, ...
        'FileFormat', 'png', ...
        'Size', 8*[1 0.75]);
    close(figH);
end



%%%%%%%%%%%%%%%%
%% overall pix AVERAGE timecourse
pMF = tcFilter(pM, ...
               'FrameTimeMs', le.FrameTimeMs, ...
               'DoHighpassFilter', true, ...
               'FilterType', 'interp', ...
               'BaselineFrames', as.interpFrNs, ...
               'FilterInterpFrames', as.interpFrNs, ...
               'DoDFOF', true, ...
               'FrameNsToAverage', as.tsAvg.frameNsToAverage, ...
               'DoPostSmooth', false);


xvals = ( (1:as.tsAvg.nFramesInBlock) - as.tsAvg.blockStimNs(1)).*(le.FrameTimeMs/1000);
figH = figure;
hold on;
plot(xvals, pMF, '-x');
vert_lines(xvals(as.tsAvg.blockStimNs));
title(sprintf('%s: %s %s', ...
              'all pix average mean', le.ExptName, le.SeriesName));
xlabel('Time (s)');
ylabel('dF/F_0');

if doExport
    set(figH, 'Visible', 'off');
    outName = fullfile(outputDir, [filePrefix 'all_pix_avg_mean']);
    exportfigPrint(figH, outName, ...
        'FileFormat', 'png', ...
        'Size', 8*[1 0.75]);
    close(figH);
end







%%%%%%%%%%%%%%%%
%% stim_pixmap
for iS=1:le.NStimsInSer
    tCurr = le.Current(iS);
    assert(length(le.Current) == le.NStimsInSer, 'listfile bug');
    tStimNs = as.tsAvg.stimNs(iS:le.NStimsInSer:end);
    
%     %% try to auto-compute smoothing and clim
%     if isempty(le.FOVSize)
%         le.FOVSize = 360;  % default to this
%     end
%     nAvgs = length(tStimNs) ./ le.NStimsInSer;
    
    
    % one for each stim current
    [crap smOutS figH] ...
        = make_stim_pixmap(le.greenStack, ...
                           'StimFrNs', tStimNs, ...
                           'NBaseFrames', as.preStimNs, ...
                           'DoPlot', true, ...
                           'CLim', pixmapCLim, ...
                           'DoMapSmooth', true, ...
                           'MapSmoothSigma', mapSmoothSigma, ...
                           'MapSmoothWidth', 3, ...
                           'ComputeWhat', 'dFOF', ...
                           'PrintStatus', false);
                           % colormap now auto-computed in m_s_p
                           % 'Colormap', cmap_posneg_yck(256, 57), ...
                       
    title(sprintf('\\bf%s\\rm: %s %s, %duA (curr %d/%d), train %sms', ...
                  'make\_stim\_pixmap', le.ExptName, le.SeriesName, ...
                  tCurr, iS, le.NStimsInSer, ...
                  le.TrainTimeMs));  % train time is a string

    cMap = get(gcf, 'Colormap');

    if doExport
        set(figH, 'Visible', 'off');
        outName = fullfile(outputDir, ...
            sprintf('%smake_stim_pixmap-curr%02d', filePrefix, iS));
        exportfigPrint(figH, outName, ...
            'FileFormat', 'png', ...
            'Size', 8*[1 0.75]);
        close(figH);
    end
end

%error('stopping before movie');
if doMovie
    %%%%%%%%%%%%%%%%
    %% PST movie
    
    avgA = stim_make_frame_ns(nFrames, le.StimEvery, 3, floor(le.StimEvery/2), ...
        le.NStimsInSer);
    % dF
    [movFr scaleRangeDF] ...
        = make_stim_avgmovie(le.greenStack, ...
        'FrameNsToAverage', avgA.frameNsToAverage, ...
        'BlockStimNs', avgA.blockStimNs, ...
        'DoPlot', false, ...
        'DFOFOutputRange', pixmapCLim, ...
        'MovieType', 'dF');
    outName = fullfile(outputDir, [filePrefix 'stim_avgmovie-dF']);
    stackWriteAvi(movFr, outName, 5, gray(256));
    
    % dFOF
    [movFr scaleRangeDFOF] ...
        = make_stim_avgmovie(le.greenStack, ...
        'FrameNsToAverage', avgA.frameNsToAverage, ...
        'BlockStimNs', avgA.blockStimNs, ...
        'DoPlot', false, ...
        'MovieType', 'dFOF', ...
        'DoSmooth', true, ...
        'SmoothWidth', 1.5, ...
        'SmoothSigma', 5, ...
        'DFOFOutputRange', [-0.15, 0.4], ...
        'DFOFBlockBaseNs', 1:5);
    
    if doExport
        outName = fullfile(outputDir, [filePrefix 'stim_avgmovie-dFOF']);
        stackWriteAvi(movFr, outName, 5, cmap_posneg_yck(256,floor(0.15/0.55*256)));
    end
end


%%%%%%%%%%%%%%%%
%% runinfo.txt
% save important data
if doExport
    outTxtName = fullfile(outputDir, [filePrefix 'runinfo.txt']);
    fid = fopen(outTxtName, 'w+');
    fprintf(fid, '\n*** %s %s\n', mfilename, datestr(now));
    fprintf(fid, 'Expt: %s %s\n', le.ExptName, le.SeriesName);
    fprintf(fid, 'nRows, nCols, nFrames: %s\nstimEvery %3d\n', ...
        mat2str(size(le.greenStack)), le.StimEvery);
    fprintf(fid, 'Listfile: Current %suA, TrainTimeMs %d, FrameTimeMs %d\n', ...
        mat2str(le.Current(:)'), le.TrainTimeMs, le.FrameTimeMs);
    if ~doMovie
        fprintf(fid, 'No movies produced\n');
    else
        fprintf(fid, 'Movie scale ranges: dF %s;   dF/F: %s\n', ...
            mat2str(scaleRangeDF(:)'), mat2str(scaleRangeDFOF(:)'));
    end

    fclose(fid);
end



                                 





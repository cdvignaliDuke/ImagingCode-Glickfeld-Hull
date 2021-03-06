function figH = stackZProfiler(stack, varargin)
%stackZProfiler (calcium): interactively plot third dim (usu. time) of stack
%
%   See below for defaults/params.
%
%   Current timecourse is in figure's UserData.timecourse field
% 
%   This now works with Matlab R2007b and uses a callback mechanism to
%   receive input rather than using the main execution thread (like
%   stack_zprofiler2).  That means that the stack is held in memory in
%   the figure's UserData structure.  I think if you share the stack
%   between windows it's ok; MATLAB uses a pointer or RCU memory to share
%   multiple read copies.  But if you want to clear that stack memory,
%   hit q to remove it from the fig's userdata struct.
%
%$Id: stackZprofiler.m 370 2008-10-28 15:32:01Z histed $

uDefs = { 'ImageToSelect', [], ...
          'Colormap', cmap_green(256), ...
          'ImageAxesProps', {}, ...
          'MaskRGBColor', [0 0 1], ...
          'MaskUseLines', false, ... % if true draw lines, else pixels
          'MaskLineColor', [1 0.5 1], ... % pinkish
          'DoImageColorbar', false, ...
          'StartRoiPolyXY', [], ...
          ...
          'DoDrawTimecourse', true, ...
          'ReducePixelMethod', [], ...
          ... % tc processing options
          'DoProcessTimecourse', true, ...
          'DoHighpassFilter', true, ...
          'FilterType', 'firls', ...
          'FilterInterpFrames', [], ...
          'HighpassHz', 0.05, ...
          'BaselineFrames', [], ...
          'DoDFOF', true, ...
          'FrameNsToAverage', [], ...
          'FrameNsToExtract', [], ...
          'DoDrawMeanOverCells', false, ...
          'FrameTimeMs', [], ...
          'DoPostSmooth', false, ...
          'SmoothMethod', 'lowess', ...
          'SmoothSpan', 5, ...
          ...  % display options
          'ZeroFrame', [], ...
          'VertLineFrameNs', [], ...
          'VertShadeFrameNs', {}, ... % cell vect, one for each shading
          'VertShadeColors', {}, ...  % cell vect, one for each shading       
          'CallFHOnTimecourse', [], ...
          'TimecourseAxesProps', {} };

%%%%%%%%%%%%%%%%
% Is this a callback?
if length(stack)==1 
    assert(ishandle(stack), 'stack error');
    h = stack;
    eventdata = varargin{1};
    idStr = varargin{2};
    
    subDealWithCallback(h, eventdata, idStr);
    return
end

%%%%%%%%%%%%%%%%
% Not a callback, called from the command line.
uo = stropt2struct(stropt_defaults(uDefs, varargin));

% check opts 
if ~strcmp(version('-release'), '2007b')
    fprintf(1, 'This version, stackZprofiler, is for MATLAB7.5 (R2007b) only');
end

%%% Setup initial frame

[nRows,nCols,nFrames] = size(stack);

% disp a warning if other active figures found
fH = findobj(0, 'Tag', 'ZprofActiveFigure');
if ~isempty(fH)
    fprintf(1, '%s: Active figures already exist; using lots of memory?\n', ...
            mfilename);
end

figH = figure;

imgAxH = subplot(1,2,1);
hold on;
set(imgAxH, 'Tag', 'ZprofImageAxes');
set(imgAxH, 'YDir', 'reverse', ...
            'XLim', [1 nCols], ...
            'YLim', [1 nRows], ...
            'OuterPosition', [0 0 0.5 1], ...
            'DataAspectRatio', [1 1 1]);
if ~isempty(uo.ImageAxesProps)
    set(imgAxH, uo.ImageAxesProps{:});
end

tcAxH = subplot(1,2,2);
hold on;
set(tcAxH, 'Tag', 'ZprofTimecourseAxes', ...
           'NextPlot', 'add', ...
           'Box', 'on');


axes(imgAxH);
if isempty(uo.ImageToSelect)
    imgToSelect = stack(:,:,1);
else
    % use the supplied image
    imgToSelect = uo.ImageToSelect;
end

% draw selection image
[nIRows nICols nIFrames] = size(imgToSelect);
assert(nIFrames == 1 || nIFrames == 3);
assert(all([nIRows nICols] == [nRows nCols]));
iH = imagesc(imgToSelect);
set(iH, 'HitTest', 'off');
set(imgAxH, 'YDir', 'reverse');

% set up selection image axes
colormap(uo.Colormap);
set(imgAxH, 'XTick', [], ...%1 nCols], ...
            'YTick', []); %[1 nRows]);
if uo.DoImageColorbar; colorbar; end

%% get initial polygon and mask: from input or by manually selecting
if isempty(uo.StartRoiPolyXY)
    tStr = 'Draw ROI, double-click when finished';
    fprintf(1, [tStr '\n']);
    title(tStr);
    [bwMask,xp,yp] = roipoly;
else
    assert(iscell(uo.StartRoiPolyXY) && prod(size(uo.StartRoiPolyXY)) == 2, ...
           'StartROIPolyXY must be a cell; first el is x, 2nd y');
    [xp,yp] = deal(uo.StartRoiPolyXY{:});
    bwMask = poly2mask(xp,yp,nRows,nCols);
end
    
% draw region
tc = subDrawRegionAndTc(bwMask, imgAxH, tcAxH, stack, uo, xp, yp);

%% define callbacks
% first image axes click
fH = get(imgAxH, 'Children');
set(fH, 'HitTest', 'off');
ud.xp = xp; 
ud.yp = yp; 
ud.nRows = nRows; 
ud.nCols = nCols;
ud.imgAxH = imgAxH; 
ud.tcAxH = tcAxH;
ud.figH = figH;
ud.stack = stack;
ud.uo = uo;
ud.timecourse = tc;
set(figH, 'UserData', ud);
set(imgAxH, 'ButtonDownFcn', {@stackZprofiler, 'Click'});
% now keypress
plotedit off;  % matlab bug: sometimes on for new figs, preventing keypressfcn
set(figH, 'KeyPressFcn', {@stackZprofiler, 'Keypress'});

% reminder to inactivate
set(figH, 'Tag', 'ZprofActiveFigure');
axes(imgAxH);
title(sprintf(['Click to move ROI.  ''q'' to quit and clear stack memory\n' ...
              '''n'' to draw a new ROI']));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bwMask, xpOut, ypOut] = subReselectRoi(imgAxH, xp, yp, nRows, nCols)

newPt = get(imgAxH, 'CurrentPoint');
figH = get(imgAxH, 'Parent');
currentChar = get(figH, 'CurrentCharacter');
newPt = newPt(1,1:2);  % reduce to xy  (comes out as 2x3 mat
if any(newPt > [nCols nRows] | newPt < [1 1])
    % error in selection, leave xp,yp as given
    keyboard
    tStr = '** Selection outside image, try again';
    fprintf(1, [tStr '\n']);
    title(tStr);
    xpOut = xp;
    ypOut = yp;
else
    % move
    roiCent = [mean(xp) mean(yp)];
    deltaC = newPt - roiCent;
    
    xpOut = xp + deltaC(1);
    ypOut = yp + deltaC(2);
end
bwMask = poly2mask(xpOut,ypOut,nRows,nCols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tcOut] = subDrawRegionAndTc(bwMask, imgAxH, tcAxH, stack, uo, ...
                                      xp, yp)

%%% check args
nPix = sum(bwMask == 1);
if nPix == 0
    tStr = '0 pixels selected; reselect';
    fprintf(1, [tStr '\n']);
    title(tStr);
    tcOut = [];
    return
end
if isempty(uo.FrameTimeMs)
    isXAxisTime = false;  % frames
else
    isXAxisTime = true; 
end

fprintf(1, 'Computing timecourse and drawing... ');
tc = stackGetTimeCourses(stack, bwMask, uo.ReducePixelMethod);
if uo.DoProcessTimecourse
    tcF = tcFilter(tc, ...
                   'FrameTimeMs', uo.FrameTimeMs, ...
                   'HighpassHz', uo.HighpassHz, ...
                   'FilterType', uo.FilterType, ...
                   'BaselineFrames', uo.BaselineFrames, ...
                   'DoHighpassFilter', uo.DoHighpassFilter, ...
                   'FilterInterpFrames', uo.FilterInterpFrames, ...
                   'DoDFOF', uo.DoDFOF, ...
                   'FrameNsToAverage', uo.FrameNsToAverage, ...
                   'DoPostSmooth', uo.DoPostSmooth, ...
                   'SmoothMethod', uo.SmoothMethod, ...
                   'SmoothSpan', uo.SmoothSpan);
else
    tcF = tc;
end

% pull out multiple timecourses if requested
if ~isempty(uo.FrameNsToExtract)
    assert(isempty(uo.FrameNsToAverage), ...
           'Cannot specify both FrameNsTo Average and Extract');
    tcE = tcF(uo.FrameNsToExtract)';
    
    tcF = tcE;  % tcF is what gets plotted
end

[nRows nCols] = size(bwMask);
[nFrames nCells] = size(tcF);

% check for correct shape of averaging matrix
if ~isempty(uo.FrameNsToExtract) || ~isempty(uo.FrameNsToAverage)
    if ~isempty(uo.ZeroFrame) ...
        && ( (uo.ZeroFrame > nFrames || uo.ZeroFrame < 1) ...
             || (any(uo.VertLineFrameNs > nFrames)))
        warning('Timecourse is size (%d,%d) (nCells,nFrames): needs transpose?', ...
            nCells, nFrames);
    end
end

       




% remove old ROI
fH = findobj(imgAxH, 'Tag', 'zprofRegionMask');
delete(fH);

% draw ROI on image plot  (always done)
axes(imgAxH);

if uo.MaskUseLines
    % draw an outline
    mH = findobj(gcf, 'Tag', 'zprofRegionLines');  % remove old
    delete(mH);
    mH = plot(xp, yp, ...
        'Color', uo.MaskLineColor, ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 3, ...
        'Tag', 'zprofRegionLines', ...
        'HitTest', 'off');
else
    % draw an image mask (pixels)
    iH = roiDrawOnImage(bwMask, uo.MaskRGBColor);
    set(iH, 'Tag', 'zprofRegionMask', ...
        'HitTest', 'off');
end


% draw timecourses if requested
axes(tcAxH);

cla;


if uo.DoDrawTimecourse
    xvals = 1:nFrames;    
    if ~isempty(uo.ZeroFrame);
        xvals = xvals - uo.ZeroFrame;
    end
    if isXAxisTime
        xvals = xvals .* uo.FrameTimeMs / 1000;
    end

    pH = plot(xvals, tcF, 'x-');
    set(pH, 'Tag', 'ZprofilerTimecourse');

    if uo.DoDrawMeanOverCells
        assert(nCells>1, 'Asked for cell mean but only one cell');

        % draw a mean
        mH = plot(xvals, mean(tcF,2)); 
        set(mH, 'LineWidth', 3, ...
                'Marker', 'x', ...
                'Color', 'k');
        set(pH, 'Marker', 'none');  % no marker for individ lines if mean 
        anystack(mH, 'top');
        
        
    end
    
    if isXAxisTime
        xlabel('Time (s)');
    else
        xlabel('Frame number');
    end
    ylabel('dF/F');
    xLim = 1.01*[min(xvals) max(xvals)];
    set(tcAxH, 'XLim', xLim);

    if ~isempty(uo.TimecourseAxesProps)
        set(tcAxH, uo.TimecourseAxesProps{:});
    end
    
    if isXAxisTime
        vert_lines(uo.VertLineFrameNs .* uo.FrameTimeMs/1000 );
    else
        vert_lines(uo.VertLineFrameNs);
    end
    
    if ~isempty(uo.VertShadeFrameNs)
        nShades = length(uo.VertShadeFrameNs);
        for iS = 1:nShades
            tShadeNs = uo.VertShadeFrameNs{iS};
            if isXAxisTime
                tShadeNs = tShadeNs .* uo.FrameTimeMs/1000;
            end
            tColor = uo.VertShadeColors{iS};
            vert_shade(tShadeNs, tColor);
        end
    else
        assert(isempty(uo.VertShadeColors), ...
               'must pass both VertShade FrameNs and Colors');
    end
    
    
end

% call function on timecourse if requested
if ~isempty(uo.CallFHOnTimecourse)
    feval(uo.CallFHOnTimecourse, tcF);
end

% output
tcOut = tcF;

% leave image axes active
axes(imgAxH);

fprintf(1, 'done.\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subDealWithCallback(h, eventdata, idStr)

switch idStr
      
  case 'Click'
    %%% mouse click on image axes.  Move ROI and update timecourses
    imgAxH = h;
    figH = get(imgAxH, 'Parent');
    ud = get(figH, 'UserData');
    % move ROI
    [bwMask, ud.xp, ud.yp] = subReselectRoi(imgAxH, ud.xp, ud.yp, ud.nRows, ud.nCols);
    % draw new one
    tc = subDrawRegionAndTc(bwMask, imgAxH, ud.tcAxH, ud.stack, ud.uo, ud.xp, ud.yp);
    
    % reload callback
    ud.timecourse = tc;
    set(figH, 'UserData', ud);

  case 'Keypress'
    %%% keypress in figure window
    figH = h;
    ud = get(figH, 'UserData');
    key = get(figH, 'CurrentCharacter');
    %key = eventdata.Key  % unreliable for some reason
    switch key
      case 'q'
        %% make window inactive
        % disable callbacks
        set(ud.imgAxH, 'ButtonDownFcn', []);
        set(figH, 'KeyPressFcn', []);
        % clear memory
        ud.stack = [];
        set(figH, 'UserData', ud);
        % remove tag
        set(figH, 'Tag', 'ZprofInactiveFigure');
        
        % display info
        axes(ud.imgAxH)
        title('');
        fprintf(1, 'Zprofiler: made inactive, stack memory cleared.\n');

      case 'n'
        %% draw a new ROI
        tStr = 'Draw new ROI, double-click when finished';
        fprintf(1, [tStr '.\n']);
        title(tStr);
        axes(ud.imgAxH);
        [bwMask,ud.xp,ud.yp] = roipoly;
        tc = subDrawRegionAndTc(bwMask, ud.imgAxH, ud.tcAxH, ud.stack, ud.uo, ud.xp, ud.yp);
        
        % resave info
        ud.timecourse = tc;
        set(gcf, 'UserData', ud);
      otherwise
        fprintf(1, 'Invalid key press: %s\n', eventdata.Key);
    end

end  % switch idStr

return



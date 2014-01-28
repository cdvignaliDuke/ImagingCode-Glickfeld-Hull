function stack_stepper(stack, varargin)
%
%$Id: stack_stepper.m 81 2008-02-23 21:50:19Z histed $

%%%%%%%%%%%%%%%%
% Is this a callback?
if prod(size(stack))==1 
    assert(ishandle(stack), 'param error');
    h = stack;
    eventdata = varargin{1};
    idStr = varargin{2};
    
    subDealWithCallback(h, eventdata, idStr);
    return
end

%%%%%%%%%%%%%%%%
% Not a callback, called from the command line, draw first frame

if nargin < 2, startFr = 1; else startFr = varargin{1}; end

[nRows,nCols,nFrames,nPlanes] = size(stack);

figH = stack_newplotfig;
axH = axes;

currFr = startFr;

subDrawFrame(currFr, stack);

%% register callbacks
% keypress
set(figH, 'KeyPressFcn', {@stack_stepper, 'KeyPressed'});

% save userdata
ud.axH = axH;
ud.figH = figH;
ud.stack = stack;
ud.currFr = currFr;
set(figH, 'UserData', ud);

return
%%% end

%%%%%%%%%%%%%%%%

function subDealWithCallback(h, eventdata, idStr)

figH = h;
ud = get(figH, 'UserData');
[nRows,nCols,nFrames,nPlanes] = size(ud.stack);

switch idStr
  case 'KeyPressed'
    switch eventdata.Key
      case 'comma'
        %% next frame
        ud.currFr = ud.currFr+1;
        if ud.currFr > nFrames, ud.currFr = 1; end % deal w overflow
        subDrawFrame(ud.currFr, ud.stack);
      case 'period'
        %% next frame
        ud.currFr = ud.currFr-1;
        if ud.currFr < 1, ud.currFr = nFrames; end % deal w overflow
        subDrawFrame(ud.currFr, ud.stack);
      case 'shift'
        % ignore
      otherwise
        fprintf(1, '%s (fig %d): unknown key %s\n', ...
                mfilename, ud.figH, eventdata.Key);
    end
end

% save UserData
set(figH, 'UserData', ud);

%%%%%%%%%%%%%%%%

function subDrawFrame(frameNum, stack);

imshow(squeeze(stack(:,:,frameNum,:)), 'InitialMagnification', 'fit');

subTitle(frameNum, stack);
colormap(gray);
set(gca, 'DataAspectRatio', [1 1 1]);
drawnow;

%%%%%%%%%%%%%%%%

function subTitle(frameNum, stack);
[nRows,nCols,nFrames,nPlanes] = size(stack);

stWD = whos('stack');

titStr = sprintf('\\bf%s\\rm: %sb: %dx%dpix, %d/%d frames', ...
                 strrep(mfilename, '_', '\_'), ...
                 num2str_metric(stWD.bytes,3), nRows, nCols, ...
                 frameNum, nFrames);
title(titStr);



function [x,y,selectionType] = getAPointbroken(axH)
%GETAPOINT (ca-tools): retrieve a single point from an axes
%
%   [x y] = getAPoint(axH);
%
%   Waits for a click OR a keystroke.  If a key is pressed,
%   [NaN key] is returned (i.e. x = NaN)
%
%   See also: GETPTS (more complicated version of this)
%$Id: getAPointSansCallbackButBroken.m 205 2008-05-02 03:37:09Z histed $

global GETAPOINT_DATA;

if nargin < 1, axH = gca; end


figH = get(axH, 'Parent');

% save state
state = uisuspend(figH);
    
% set up figure
set(figH, 'Pointer', 'crosshair');

%% wait for click
try
    wasKey = waitforbuttonpress;
catch
    if ishandle(figH)  % not closed
        uirestore(state);
    end
    rethrow(lasterror);
end
    
%% got a click, return it  
if wasKey
    key = get(figH, 'CurrentCharacter');    
    x = NaN;
    y = key;
    selectionType = NaN;
    return
else
    % click
    tPt = get(axH, 'CurrentPoint');
    x = tPt(1,1);
    y = tPt(1,2);
    selectionType = get(figH, 'SelectionType');
end











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% not used now, using 'cross' pointer
function [pointerShape, pointerHotSpot] = subCreatePointer

pointerHotSpot = [8 8];
pointerShape = [ ...
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
              1   1   1   1   1   1   1   2   1   1   1   1   1   1   1   1
              2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2
              1   1   1   1   1   1   1   2   1   1   1   1   1   1   1   1
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN   1   2   1 NaN NaN NaN NaN NaN NaN NaN
            NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];


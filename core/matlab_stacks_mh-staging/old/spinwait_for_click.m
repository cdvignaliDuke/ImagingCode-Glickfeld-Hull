function [pt, character] = spinwait_for_click(objH, doWaitForCharAlso)
%
%$Id: spinwait_for_click.m 97 2008-03-17 20:14:13Z histed $

if nargin < 1, objH = gca; end
if nargin < 2, doWaitForCharAlso = true; end  % otherwise just click

figH = gcf;

drawnow;
ptOnEntry = get(objH, 'CurrentPoint');
charOnEntry = get(figH, 'CurrentCharacter');


while 1
    drawnow;
    currPt = get(objH, 'CurrentPoint');
    currChar = get(figH, 'CurrentCharacter');
    
    if any(currPt(:) ~= ptOnEntry(:))
        pt = currPt;
        character = []; 
        return
    elseif doWaitForCharAlso && any(currChar ~= charOnEntry)
        pt = [];
        character = currChar;
        return
    end  % else keep spinning
end

    

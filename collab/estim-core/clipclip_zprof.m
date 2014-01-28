function clipclip_zprof(figH, tcAxSize, tcFormat)
%
%
% copy current figure but use screen snapshot
%
%$Id: clipclip_zprof.m 95 2008-03-16 03:47:36Z histed $

if nargin < 1 || isempty(figH), figH = gcf; end
if nargin < 2 || isempty(tcAxSize), tcAxSize = 6*[1 1]; end
if nargin < 3, tcFormat = 'pdf'; end

axH = findobj(figH, 'Type', 'axes');

ud = get(gcf, 'UserData');

% hide title: not req because bounding box is tight around axes
%titH = get(ud.imgAxH, 'Title');
%set(titH, 'Visible', 'off');
clipclipimage(ud.imgAxH);
%set(titH, 'Visible', 'on');


% dup axes
newFigH = copyobj(figH, 0);
delete(get(newFigH, 'Children'));
newAxH = copyobj(ud.tcAxH, newFigH);
axSize = get(newAxH, 'OuterPosition');
set(newAxH, 'OuterPosition', [0 0 1 1]);

clipclip(newFigH, tcAxSize, tcFormat);
delete(newFigH);
return

% $$$ set(figH, 'Color', 'w');
% $$$ 
% $$$ 
% $$$ set(axH, 'FontSize', 14);
% $$$ set(0, 'ShowHiddenHandles', 'on');
% $$$ tH = findobj(figH, 'Type', 'text')
% $$$ set(tH, 'FontSize', 16);
% $$$ set(0, 'ShowHiddenHandles', 'off');
% $$$ 
% $$$ fprintf(1, 'Now press Cmd-shift-4 SPC then click on the window\n');
% $$$ 
% $$$ %unix(sprintf('convert %s %s', [fullF '.png'], [fullF '.wmf']));

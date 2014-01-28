function [figH] = make_roi_with_pixelcount(img, xpIn, ypIn)
%
%  Tells you how many pixels are in ROI while drawing.
%
%  To get out data:
%  ud = get(gcf, 'UserData');
%  ud.bwMask, ud.xp, ud.yp contain the logical mask, x points and y
%  points, successively.
%
%  To move polygon vertices from command line:
%  rH = findobj(gcf, 'Tag', 'RoiPixelcount');
%  pos = get(rH);
%  change pos and then 
%  set(pH, pos);  
%  pixelcount will automatically update
% 
%$Id: make_roi_with_pixelcount.m 79 2008-01-17 17:05:40Z histed $

[nRows nCols nPlanes] = size(img);
assert(nPlanes == 3 || nPlanes == 1);
if nargin < 2, xpIn = floor(nCols*[1 1 2 2]/3); end
if nargin < 3, ypIn = floor(nRows*[1 2 2 1]/3); end

%%%

figH = figure;

drawnow;
set(gcf, 'Selected', 'off', ...
         'ToolBar', 'figure');
plotedit('off');  % not sure why this is getting turned on


imH = imagesc(img);
axH = get(imH, 'Parent');
set(axH, 'DataAspectRatio', [1 1 1]);

% create impoly
pH = impoly(axH, cat(2, xpIn(:), ypIn(:)));
set(pH, 'tag', 'PixelcountRoi');

%% create text area
tH = uicontrol('Parent', figH, ...
               'Style', 'text', ...
               'Units', 'normalized', ...
               'Position', [0.01 0.90 0.1 0.08], ...
               'FontSize', 8, ...
               'FontName', 'Helvetica', ...
               'HorizontalAlignment', 'left', ...
               'Tag', 'NPixText');

               
               
%% t
% register callbacks
api = iptgetapi(pH);
api.addNewPositionCallback(@subNewPosCB);

% update
subUpdateText(figH, xpIn(:), ypIn(:), nRows, nCols);

%%%%%%%%%%%%%%%%

function subNewPosCB(pos)
xp = pos(:,1);
yp = pos(:,2);

% get image axis
imH = findobj(gca, 'Type', 'Image');
if length(imH) ~= 1 
    fprintf(1, 'Error: multiple images in axes?\n');
end

[nRows nCols nPlanes] = size(get(imH, 'CData'));

subUpdateText(gcf, xp, yp, nRows, nCols);

%%%%%%%%%%%%%%%%

function subUpdateText(figH, xp, yp, nRows, nCols)

bwMask = poly2mask(xp, yp, nRows, nCols);

nPix = sum(bwMask(:));

% change text
tH = findobj(figH, 'Tag', 'NPixText');
set(tH, 'String', ...
        sprintf('nPix in area:\n%d', nPix));

ud.bwMask = bwMask;
ud.xp = xp;
ud.yp = yp;
set(figH, 'UserData', ud);


    

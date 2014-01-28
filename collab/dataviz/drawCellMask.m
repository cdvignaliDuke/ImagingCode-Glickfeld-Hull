function [figH,tH] = drawCellMask(labelMask, fov, plotWhat)
%DRAWCELLMASK
%  FIGH = DRAWCELLMASK(LABELMASK, FOV, PLOTWHAT)
%
%  plotWhat:  how to plot cells: 'outlines', 'numbers', 'patches'
%     can specify a cell array with any set of these listed i.e. 
%     {'numbers','outlines'}
%
%  if fov is specified, plot it as the underlying image.  Else use the
%  label mask itself
%
%  set showFlag to false to make figure invisible
%
%  MH 080501:  import to calcium; add numbers/outlines options
%  MH 081020:  allow multiple plotWhat options
%  VB 080825:  added showFlag
%
%$Id: figCellsLabeled.m 398 2008-11-24 17:30:12Z vincent $

if nargin < 2, fov = []; end
if nargin < 3, plotWhat = 'outlines'; end

doFovBackground = ~isempty(fov);

%%
[nRows nCols] = size(labelMask);
rp = regionprops(labelMask);

figH = gcf;

if doFovBackground
    if size(fov,3) > 1
        image(fov);
    else
        imagesc(fov);
        colormap(cmap_green);
    end
    textColor = [1 1 1];
else
%     imagesc(labelMask);
     textColor = [0 0 0];
end;

hold on;

% parse plotWhat
doNumbers = false;
doOutlines = false;
doPatches = false;
if ischar(plotWhat), plotWhat = {plotWhat}; end
if any(ismember(plotWhat, 'numbers'))
    doNumbers = true;
end
if any(ismember(plotWhat, 'outlines'))
    doOutlines = true;
end
if any(ismember(plotWhat, 'patches'))
    doPatches = true;
end
if doNumbers == false && doOutlines == false && doPatches == false
    error('No matching arguments to plotWhat specified (must be lowercase)');
end

%% draw on image

if doPatches
    hotPink = [237 40 145] / 255;  % default color
    
    bwMask = labelMask > 0;
    bounds = bwboundaries(bwMask, 8, 'noholes');
    nB = length(bounds);
    for iB = 1:nB
        patch(bounds{iB}(:,2), bounds{iB}(:,1), hotPink, ...
             'Tag', 'CellPatches', ...
              'EdgeColor', 'none');
        
    end
end

if doOutlines
%     bwMask = labelMask > 0;
%     bounds = bwboundaries(bwMask, 8, 'noholes');
%     nB = length(bounds);
    
    for iB = unique(labelMask(:))'
        bwMask = labelMask == iB;
        bound = bwboundaries(bwMask,8,'noholes');
        
        plot(bound{1}(:,2), bound{1}(:,1), 'w', ...
             'Tag', 'CellOutlines');
        
%         plot(bounds{iB}(:,2), bounds{iB}(:,1), 'w', ...
%              'Tag', 'CellOutlines');
    end
end

tH=[];
if doNumbers
    nCells = max(labelMask(:));
    labList = 1:nCells;
    for iL=labList
        tH(iL) = text(rp(iL).Centroid(1), rp(iL).Centroid(2), ...
                      sprintf('%d',iL), ...
                      'HorizontalAlignment', 'center', ...
                      'VerticalAlignment', 'middle', ...
                      'Tag', 'CellNumbers', ...
                      'Color', textColor);
    end
end

if ~doFovBackground
    % make white background
    cMap = get(gcf, 'Colormap');
    cMap(1,:) = [1 1 1];
    set(gcf, 'Colormap', cMap);
end

set(gca, 'TickDir', 'out', 'xtick',[],'ytick',[],...
         'Box', 'on');
     
%          'XTick', [1 nCols], ...
%          'YTick', [1 nRows], ...

set(gca,'plotboxaspectratio',[nCols nRows 1]);

% export
%clipclip(gcf, 6*[1 1]);

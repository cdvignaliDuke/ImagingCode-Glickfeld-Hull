function lastShift = find_phase_adjust_fromimage(tImage, cMap)
%
% example call:
%    find_phase_adjust(fov)
%
%$Id: find_phase_adjust_fromimage.m 78 2008-01-10 15:38:51Z histed $

if nargin < 2, cMap = []; end

[nRows nCols nPlanes] = size(tImage);

figH = figure;
axH = axes('YDir', 'reverse', ...
           'XLim', [1 nCols], ...
           'YLim', [1 nRows], ...
           'DataAspectRatio', [1 1 1]);
hold on;

nShift = 0;
while 1
    shiftImage = adjust_phase(tImage, nShift);

    cla;
    iH = imagesc(shiftImage);    
    if ~isempty(cMap)
        colormap(cMap);
    end
    
    title(make_plot_title(sprintf(': phase shift %d', nShift)));    

    drawnow;

    lastShift = nShift;
    nShift = input('next nShift: ');
    if isempty(nShift)
        fprintf(1, 'Done, last shift %d\n', lastShift);
        return
    end
end







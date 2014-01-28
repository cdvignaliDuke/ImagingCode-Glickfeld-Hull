function figH = make_stim_cellmaskfig(greenStack, labNeur, varargin)
%
%   function figH = make_stim_cellmaskfig(greenStack, labNeur, STROPTS)
%
%$Id$

defs = { 'StimFrNs', [], ...
         'StimFrPeakOffset', 0, ...
         'BasePreOffset', 0, ...
         'PostFrameNs', 1, ...
         'NBaseFrames', 10, ...
         'CLim', [0 0.20], ...
         'Colormap', hot, ... 
         'DoColorbar', true, ...
         'MinDenomF', 0, ...
         'MinDenomF', 40 };

     
uo = stropt2struct(stropt_defaults(defs, varargin));

%%%

[nRows nCols nFrames] = size(greenStack);

% get pixmap
outS = make_stim_pixmap(greenStack, ...
    'StimFrNs', uo.stimFrNs, ...
    'StimFrPeakOffset', uo.StimFrPeakOffset, ...
    'BasePreOffset', uo.BasePreOffset, ...
    'PostFrameNs', uo.PostFrameNs, ...
    'NBaseFrames', uo.NBaseFrames, ...
    'DoPlot', false, ...
    'DoMapSmooth', false, ...
    'ComputeWhat', 'dFOF', ...
    'MinDenomF', uo.MinDenomF, ...
    'PrintStatus', false);

% compute cell means
nCells = max(labNeur(:));
outImg = outS*0;
for iC=1:nCells
    tVs = outS(labNeur==iC);
    tM = mean(tVs(:));
    outImg(labNeur==iC) = tM;
end

%% make plot
figH = figure;
imagesc(outImg);
set(gca, 'CLim', uo.CLim);
colormap(uo.Colormap);
axis square
colorbar;

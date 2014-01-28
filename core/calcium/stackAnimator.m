function stackAnimator(stack, varargin)
%STACKANIMATOR Show stack as a movie in a matlab fig
% STACKANIMATOR(STACK, VARARGIN)
%   if 4d, the fourth dimension is the RGB dim
%
%   Opts:
%   'FramesPerSec': scalar
%   'DoColobar': boolean
%   'edit stackAnimator' to see full list of options.  
%  
%   See also: STACKZPROFILER
%
%$Id: stackAnimator.m 500 2009-04-21 20:42:16Z vincent $

defs = { 'FramesPerSec', 15, ...
         'Colormap', gray(256), ...
         'CLim', [], ...
         'DoColorbar', false, ...
         'AxesProps', {}, ...
         'FigProps', {} };

uo = stropt2struct(stropt_defaults(defs, varargin));

[nRows,nCols,nFrames,nChans] = size(stack);
tic;

figH = figure;
if ~isempty(uo.FigProps)
    set(figH, uo.FigProps{:});
end
axH = axes;

%% sample pix to find CLim
if isempty(uo.CLim)
    t = toc;
    sampPix = stack(floor(linspace(1,nRows,10)), floor(linspace(1,nCols,10)),:);
    sampPix = sampPix(:);
    allMin = min(sampPix);
    allMax = max(sampPix);
    cLim = [allMin allMax];
    if (toc-t) > 3
        disp(toc-t);
        warning('That took too long: fix?');
    end
else
    cLim = uo.CLim;
end



% loop by default
currFr = 1;
iH = [];
fprintf(1, '%s: hit q to exit ... ', mfilename);

while true
    
    t2 = toc;
    if ~isempty(iH), delete(iH); end
    axes(axH);    
    iH = imagesc(squeeze(stack(:,:,currFr,:)));
    %set(iH, 'Parent', axH);
    set(axH, 'OuterPosition', [0 0 1 1], ...
             'DataAspectRatio', [1 1 1], ...
             'CLim', cLim);
    subTitle(currFr, stack, axH);
    set(figH, 'Colormap', uo.Colormap);
    set(axH, 'DataAspectRatio', [1 1 1], ...
             uo.AxesProps{:});
    if uo.DoColorbar
        colorbar('peer', axH);
    end
    drawnow;
    
    currFr = currFr+1;
    if currFr > nFrames
        currFr = 1;
    end

% $$$     % check for keypress
% $$$     keyPressed = ~isempty(get(gcf, 'CurrentCharacter'));
% $$$     mouseClicked = ~all(get(gcf, 'CurrentPoint') == [0 0]);
% $$$     if keyPressed || mouseClicked
% $$$         fprintf('%s: Key pressed/mouse clicked, stopping animation\n', ...
% $$$                 mfilename);
% $$$         break
% $$$     end
    
    % wait until next fr start
    elS = toc - t2;
    pause(max(1/uo.FramesPerSec - elS, 0));  % subtr execution time
    
    if lower(get(figH, 'CurrentCharacter')) == 'q';
        fprintf(1, 'done\n');
        break;
    end
end


%%%%%%%%%%%%%%%%

function subTitle(frameNum, stack, axH);
[nRows nCols nFrames nColors] = size(stack);

stWD = whos('stack');

titStr = sprintf('\\bf%s\\rm: %sb: %dx%dpix, %d / %d frames', ...
                 strrep(mfilename, '_', '\_'), ...
                 num2str_metric(stWD.bytes,3), nRows, nCols, ...
                 frameNum, nFrames);
title(axH, titStr);


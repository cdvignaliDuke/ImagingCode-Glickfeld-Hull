function outStack = stack_timestamper(stack, frameHz, unitStr)
%
%$Id: stack_timestamper.m 58 2007-11-15 21:38:14Z histed $

[nRows,nCols,nFrames,nPlanes] = size(stack);

if nargin < 3, unitStr = 'seconds'; end

pixStart = [nCols-10 10];
textProps = { 'Color', 'w', ...
              'HorizontalAlignment', 'right', ...
              'VerticalAlignment', 'bottom', ...
              'FontName', 'FixedWidth', ...
              'FontSize', 16, ...
              'FontWeight', 'bold' };

textStrs = {}
for iF = 1:nFrames
    frameTimeS = (iF-1) ./ frameHz;

    switch unitStr
      case 'seconds'
        textStrs{iF} = sprintf('%6.1f sec', roundto(frameTimeS,1));
      case 'minutes'
        %textStrs{iF} = sprintf('%6.1f min', roundto(frameTimeS/60,1));
        textStrs{iF} = sprintf('%6g min', chop(frameTimeS/60,2));
      otherwise
        error('not implemented yet');
    end
end

outStack = stack_annotator(stack, textStrs, 1:nFrames, pixStart, textProps);
    
    

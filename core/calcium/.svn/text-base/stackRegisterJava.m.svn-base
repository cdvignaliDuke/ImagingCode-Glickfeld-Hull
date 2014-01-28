function regStack = stackRegisterJava(stack,target)
%STACKREGISTERJAVA (calcium): register a stack (single-threaded)
% REGSTACK = STACKREGISTERJAVA(STACK)
% REGSTACK = STACKREGISTERJAVA(STACK,TARGET)
%   Essentially a matlab wrapper around IJAlign_AK.java
%
%   histed 080505: create from RegMcoreSlave.m (aaron)
%   bonin 100107 added optional target argument
%$Id$

[nRows nCols nFrames] = size(stack);

if nargin < 2
    nAvg = min(100, floor(nFrames/15));
    avgFrNs = 1:nAvg;
    target = mean(stack(:,:,avgFrNs), 3);
end

fprintf(1, 'Moving source data to Java... ');
targetImg = ijarray2plus(target, 'single');
sourceStack = ijarray2plus(stack, 'single');
al = IJAlign_AK;
clear('stack');
fprintf(1, 'done.\n');

%% set up options
cropping = sprintf('%d %d %d %d', [0 0 nCols-1 nRows-1]);
landmarks=num2str(round([nCols nRows nCols nRows]/2));
%landmarks = [center,' ',center,' ',center,' ',center];
%landmarks = center;
cmdstr = sprintf(['-align -window s %s -window t %s ' ...
                  '-translation %s -hideOutput'], ...
                 cropping, cropping, landmarks);

% run it
tic;
fprintf(1, 'Doing alignment of %d frames... ', nFrames);
resultStack = al.doAlign(cmdstr, sourceStack, targetImg);
elT = toc;
fprintf(1, 'done in %ds, %4.2ffps\n', ...
        chop(elT,2), nFrames./elT);


%% decode output

% dummy object to disable automatic scaling in StackConverter
fprintf(1, 'Moving result data back to MATLAB... ');
dummy = ij.process.ImageConverter(resultStack);
dummy.setDoScaling(0); % this static property is used by StackConverter

converter = ij.process.StackConverter(resultStack); % don't use 
                                                   %ImageConverter for stacks!
converter.convertToGray8;

regStack = ijplus2array(resultStack);
fprintf(1, 'done.\n');



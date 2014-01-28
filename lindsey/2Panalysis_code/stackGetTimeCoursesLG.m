function tc = stackGetTimeCourses(stack, labelMask, reducePixelMethod)
%STACKGETTIMECOURSES Calculate pixel intensities inside multiple ROIs
%  TC = STACK_GET_TCOURSES(STACK, LABELMASK)
%  TC = STACKGETTIMECOURSES(STACK, LABELMASK, REDUCEPIXELMETHOD)
%
%   stack is an image/stack, of size (nRows,nCols,nFrames)
%       nFrames is >= 1
%   labelMask is a double or logical matrix, size (nRows,nCols)
%   
%   tc is a matrix of size (nFrames, nCells), each cell is a mean over a
%        region specified in labelMask
%     
%$Id: stackGetTimeCourses.m 398 2008-11-24 17:30:12Z vincent $

if nargin < 3 || isempty(reducePixelMethod), reducePixelMethod = 'mean'; end

[nRows nCols nFrames]=size(stack);

if nargin < 2
    labelMask = ones(nRows,nCols);
end

if any([nRows nCols] ~= size(labelMask)) 
    error('mismatched sizes of stack and labelMask');
end
labelMask = double(labelMask);  % in case logical

nCells = max(labelMask(:));

% $$$ % let the JIT compiler make this fast    
% $$$ tic;
% $$$ tc = repmat(NaN, [nCells, nFrames]);
% $$$ for iF=1:nFrames
% $$$     for iC=1:nCells
% $$$         stFr = stack(:,:,iF);
% $$$         stVec = stFr(labelMask==iC);
% $$$         tc(iC,iF) = mean(stVec(:)); 
% $$$     end
% $$$ end
% $$$ toc

% do reshape, takes more mem, but is faster
tic;
stackCols = reshape(stack, [nRows*nCols, nFrames]);
labelCols = reshape(labelMask, [nRows*nCols, 1]);

% check mask
unqLabels = setdiff(unique(labelMask(:)), 0);  % remove 0 label
if length(unqLabels) ~= max(unqLabels)
    error('Mask labels assumed contig. increasing but N=%d, max=%d', ...
          length(unqLabels), max(unqLabels));
end
% do reduce
switch reducePixelMethod
  case 'mean'
    tc = repmat(NaN, [nFrames,nCells]);
    for iC = 1:nCells
        tc(:,iC) = mean(stackCols(labelCols==iC,:),1);
    end
  case 'min'
    tc = repmat(NaN, [nFrames,nCells]);
    for iC = 1:nCells
        tc(:,iC) = min(stackCols(labelCols==iC,:));
    end
  case 'median'
    tc = repmat(NaN, [nFrames,nCells]);
    for iC = 1:nCells
        tc(:,iC) = median(stackCols(labelCols==iC,:));
    end
  case 'none'
    assert(nCells == 1, 'none method - only one cell allowed');
    nPix = sum(labelCols==1);
    tc = stackCols(labelCols==1, :)';
  otherwise
    error('invalid reduce method: %s');
end

t=toc;
if t>5
    fprintf(1, '%s: elapsed time is %g seconds.\n', mfilename, t);
end

% ************ PREVIOUS VERSION OF THE CODE *******************
%
% function tc = stackGetTimeCourses(stack,labelimage)
% %STACKGETTIMECOURSES Get timecourses from cell ROI
% % tc = stackGetTimeCourses(stack,labelimage)
% %
% % input: stack:   3D array (series of images)
% %   labelimage: 2D array ('image') with regions defined by indices: 1, 2...
% % output:  tc   timecourses.  (tDim,nLbl)
% % From get_tcoursesSY  only 1 change: changed temporary variable flatstack
% %   to stack.  In either case, it appears as if stack is passed by value
% %   and a new temporary variable is NOT made.  Also returns transpose of
% %   the original, so that time variable is fastest moving
% % tic
% % based on get_tcoursesCR by Clay Reid
% 
% [xDim,yDim,tDim]=size(stack);
% 
% nLbl=max(max(labelimage));
% 
% 
% tc=NaN(nLbl,tDim); %changed to NaN from zeros YC 3/6/06
% 
% stack=reshape(stack,xDim*yDim,tDim);
% 
% 
% for j=1:nLbl
%     cind=find(labelimage==j);
%     if ~isempty(cind) %added YC 2/15/06 to deal with relabeled masks   
%         tc(j,:)=sum(stack(cind,:))/length(cind);
%     end
% end
% tc = tc';
% % toc

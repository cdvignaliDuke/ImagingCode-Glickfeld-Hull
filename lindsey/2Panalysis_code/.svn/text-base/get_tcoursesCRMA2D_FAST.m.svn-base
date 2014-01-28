function tc = get_tcoursesCRMA2D_FAST(stack,labelimg)
%much FASTER than get_tcoursesCRMA2D.m only when mask is distributed (e.g.
% all neuropil)
% get timecourses of stack regions defined by index 3-D array: labelvolume
% input: stack:   4D array (series of volumes)
%   labelvolume: 3D array ('image') with regions defined by indices: 1, 2...
% output:  tc   timecourses.  (tDim,nLbl)
% From get_tcoursesSY  only 1 change: changed temporary variable flatstack
%   to stack.  In either case, it appears as if stack is passed by value
%   and a new temporary variable is NOT made.  Also returns transpose of
%   the original, so that time variable is fastest moving
% tic
[xDim,yDim,tDim]=size(stack);

nLbl=max(max(labelimg));


tc=NaN(nLbl,tDim); %changed to NaN from zeros YC 3/6/06

%stack=reshape(stack,xDim*yDim,tDim);


for j=1:nLbl
%    cind=find(labelimg==j);
   mask = labelimg==j;
   N = sum(sum(mask));
   if N>0 %added YC 2/15/06 to deal with relabeled masks   
       for k = 1:size(stack,3)
           tc(j,k)=sum(sum(squeeze(stack(:,:,k)).*mask))./N;
       end
   end
end
tc = tc';
% toc

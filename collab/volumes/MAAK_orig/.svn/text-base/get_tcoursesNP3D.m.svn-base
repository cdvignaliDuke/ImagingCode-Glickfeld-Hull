function tc = get_tcoursesNP3D(stack,npmasks)
% get timecourses of stack regions defined by index 3-D array: labelvolume
% input: stack:   4D array (series of volumes)
%   labelvolume: 3D array ('image') with regions defined by indices: 1, 2...
% output:  tc   timecourses.  (tDim,nLbl)
% From get_tcoursesSY  only 1 change: changed temporary variable flatstack
%   to stack.  In either case, it appears as if stack is passed by value
%   and a new temporary variable is NOT made.  Also returns transpose of
%   the original, so that time variable is fastest moving
% tic
[xDim,yDim,zDim,tDim]=size(stack);

nLbl=size(npmasks,4);


tc=NaN(nLbl,tDim); %changed to NaN from zeros YC 3/6/06

stack=reshape(stack,xDim*yDim*zDim,tDim);

for j=1:nLbl
    cind=find(npmasks(:,:,:,j)==1);
    if ~isempty(cind) %added YC 2/15/06 to deal with relabeled masks   
        tc(j,:)=sum(stack(cind,:))/length(cind);
    end
end
tc = tc';
% toc
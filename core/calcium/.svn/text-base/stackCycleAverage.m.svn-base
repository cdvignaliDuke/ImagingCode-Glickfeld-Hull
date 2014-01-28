function [av,ncycles] = stackCycleAverage(array,nframes,c)
%STACKCYCLEAVERAGE TRIAL AVERAGE
%   [AV,NCYCLES] = STACKCYCLEAVERAGE(ARRAY,NFRAMES) 
%   [AV,NCYCLES] = STACKCYCLEAVERAGE(ARRAY,NFRAMES,C) to specify
%   the class of AV

if nargin < 3
    c = class(array);
end

[ny,nx,nsamples]=size(array);
ncycles = floor(nsamples/nframes);

s = zeros(ny,nx,nframes);
base = [0:ncycles-1]*nframes;
tic;
for iFrame = 1:nframes
    s(:,:,iFrame)=mean(array(:,:,base+iFrame),3);
end
fprintf(1,'%2.1f fps\n',nframes/toc);

av = cast(s / ncycles,c);

return
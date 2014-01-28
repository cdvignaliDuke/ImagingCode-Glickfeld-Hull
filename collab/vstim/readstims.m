function frames = readstims(pathname,nStim)
%READSTIMS
%FRAMES = READSTIMS(PATHNAME,NSTIM)

% load first frame to get resolution
fn = 'stim0 frame0.bmp';

thisFrame = imread(fullfile(pathname,fn));
[ny,nx,nz]=size(thisFrame);

list = dir(fullfile(pathname,'*.bmp'))

nFrames = length(list);

nPicsPerStim = nFrames / nStim;

frames = zeros([ny,nx,nFrames],'uint8');

index = 1;
for iStim = 1:nStim
    fprintf(1,'Loading stimulus %i/%i ',iStim,nStim);
    for iPic = 1:nPicsPerStim
        fn = sprintf('stim%i frame%i.bmp',iStim-1,iPic-1);
        thisFrame  = imread(fullfile(pathname,fn));
        frames(:,:,index) = thisFrame(:,:,1);
        index = index+1;
        if mod(iPic,1000)==0
            fprintf('%i ',iPic);
        end
    end
    fprintf('\n');
end

return
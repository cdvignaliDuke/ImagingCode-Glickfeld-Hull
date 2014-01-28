function b = writestims(imgs,pathname,stims);
% b = writestims(imgs,pathname,stims);

fprintf(1,'Writing bmp files...\n');

tic;

if nargin < 3
    if iscell(imgs) 
        stims = 0:length(imgs)-1;
    end
end

if ~exist(pathname)
    mkdir(pathname);
else
    delete(fullfile(pathname,'*.bmp'));
end

total = 0;
for istim = 1:length(stims)
    
    [ny,nx,nframes]=size(imgs{istim});

    %% generate bmp files (bmp files are faster to load)
    myimage = zeros([size(imgs{istim}(:,:,1)) 3],'uint8');
    for iframe= 1:nframes
        myimage(:,:,1) = imgs{istim}(:,:,iframe);
        myimage(:,:,2) = imgs{istim}(:,:,iframe);
        myimage(:,:,3) = imgs{istim}(:,:,iframe);
        fn = sprintf('stim%i frame%i.bmp',stims(istim),iframe-1);
        pn = fullfile(pathname,fn);
        imwrite(myimage,pn,'bmp');    
    end
    total = total + nframes;
    
end % istim
    
t = toc;

s = whos('imgs');

fprintf(1,'Wrote %i bmp files to directory ''%s'' (%2.1f MB/s)\n',total, pathname,s.bytes/1024/1024/t);

return;
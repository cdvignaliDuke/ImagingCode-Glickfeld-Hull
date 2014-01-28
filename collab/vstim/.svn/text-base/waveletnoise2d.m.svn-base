function imgs = waveletnoise2d(prot, s)
%WAVELETNOISE 
% IMGS = WAVELETNOISE2D(PROT, S)

dt = s.DecimationRatio/s.RefreshRate;

total = 0;

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

% each protocol consists of multiple stimuli
for iStim = 1:prot.nstim  
	sel = find([prot.pars.n]==iStim);
    
    % each stimulus consists of multiple segments
    for iSeg= 1:length(sel);
        fprintf(1,'Stimulus %i ',iStim);
        pars = prot.pars(sel(iSeg));
        J = log2(pars.npix)-1; % number of scales (excl lowest freq)
        
        rand('twister',pars.seed);
        
        inverse = zeros([pars.npix*2,pars.npix*2,pars.len]);
        
        for iFrame = 1:pars.len
            
            if mod(iFrame,100)==0
               fprintf('.'); 
            end
            
            % generate 4x the area and crop to avoid wrapping
            coefs = dualtree2dcoefs([pars.npix*2,pars.npix*2],J,pars.p0,pars.k,[pars.Jmin pars.Jmax]);
            inverse(:,:,iFrame)  = idualtree2D(coefs, J, Fsf, sf);
        end
        fprintf('\n');
        
        y = inverse(1:pars.npix,1:pars.npix,:);
%         lims = prctile(y(:),[.5 99.5]);

        sigma = std(y(:));
        lims = pars.range*[-sigma sigma];

        imgs{iStim} = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);
        imgs{iStim}(:,:,end) = 128;
    end
end

return;

%%
dirs = frGetDirs;
addpath(fullfile(dirs.svn,'third/dualtree'));

s = screen;
prot = readprot('prots/waveletnoise2d.prot');
imgs = waveletnoise2d(prot,s);
% writetiff(imgs{1},'test.tif');

%%
imagesc(imgs{1}(:,:,1));axis square; colormap gray;


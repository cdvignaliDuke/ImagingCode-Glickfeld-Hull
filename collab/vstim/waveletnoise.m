function imgs = waveletnoise(prot, s)
%WAVELETNOISE 
% IMGS = WAVELETNOISE(PROT, S)
% 
% DEPRECATED. USE waveletnoise3d instead!!!'

warning('DEPRECATED. USE waveletnoise3d instead!!!');
warning('DOES NOT INITIALIZE SEED!!!');

dt = s.DecimationRatio/s.RefreshRate;

total = 0;

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

% each protocol consists of multiple stimuli
for iStim = 1:prot.nstim  
	sel = find([prot.pars.n]==iStim);
    
    % each stimulus consists of multiple segments
    for iSeg= 1:length(sel);
        fprintf(1,'Stimulus %i\n',iStim);
        pars = prot.pars(sel(iSeg));
        J = log2(pars.npix)-1; % number of scales (excl lowest freq)
        % generate 4x the area and crop to avoid wrapping
        coefs = dualtree3dcoefs([pars.npix*2,pars.npix*2,pars.len],J,pars.seed,pars.p0,pars.k,[pars.Jmin pars.Jmax]);
        inverse  = idualtree3D(coefs, J, Fsf, sf);
        cropped = inverse(1:pars.npix,1:pars.npix,:);
        lims = prctile(cropped(:),[.5 99.5]);
        imgs{iStim} = uint8((cropped*pars.c-lims(1))/(lims(2)-lims(1))*255);   
        imgs{iStim}(:,:,end) = 128;
    end
end

return;

%%
dirs = frGetDirs;
addpath(fullfile(dirs.svn,'third/dualtree'));

s = screen;
prot = readprot('prots/waveletnoise.prot');
imgs = waveletnoise(prot,s);


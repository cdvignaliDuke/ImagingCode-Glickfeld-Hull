function imgs = waveletnoise3d(prot, s)
%WAVELETNOISE3D 
% IMGS = WAVELETNOISE3D(PROT, S)

% dt = s.DecimationRatio/s.RefreshRate;
% 
% This is the stimulus used for the paper entitled:
% Robust characterization of receptive fields in mouse visual cortex

% each protocol consists of multiple stimuli
for iStim = 1:prot.nstim  
	sel = find([prot.pars.n]==iStim);
    
    % each stimulus consists of multiple segments
    for iSeg= 1:length(sel);
        fprintf(1,'Stimulus %i\n',iStim);
        pars = prot.pars(sel(iSeg));

        imgs{iStim} = waveletrnd3(pars);
        imgs{iStim}(:,:,end) = 128;
    end
end

return;

%%
dirs = frGetDirs;
addpath(fullfile(dirs.svn,'third/dualtree'));

s = screen;
prot = readprot('prots/waveletnoise.prot');
imgs = waveletnoise3d(prot,s);


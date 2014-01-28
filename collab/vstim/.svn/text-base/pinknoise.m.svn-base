function imgs = pinknoise(prot, s)
%PINKNOISE Bandpass noise as described in Niell and Stryker 2008
% IMGS = PINKNOISE(PROT, S)
%
% PROT is struct describing stimuli to be generated
% S is SCREEN structure
%
% see READPROT, SCREEN
% see PINKRND, PINKNOISE2, PINKRND2

% change log
%
% 09/07/23 vb move part of code to pinkrnd
% 

refreshRate = s.RefreshRate/s.DecimationRatio;

% each protocol consists of multiple stimuli
for iStim = 1:prot.nstim  
	sel = find([prot.pars.n]==iStim);
    
    % each stimulus consists of multiple segments
    for iSeg= 1:length(sel);
        fprintf(1,'Stimulus %i\n',iStim);
        pars = prot.pars(sel(iSeg));
        
        r = pinkrnd(pars,refreshRate);

        imgs{iStim} = r;   
        imgs{iStim}(:,:,end) = 128;
    end
end

return;

%%
dirs = frGetDirs;

s = screen;
prot = readprot('prots/pinknoise.prot');
imgs = pinknoise(prot,s);
writetiff(imgs{1},'test.tif');


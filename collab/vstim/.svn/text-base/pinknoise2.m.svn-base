function imgs = pinknoise2(prot, s)
%PINKNOISE2 Bandpass noise stimulus (second generation)
% IMGS = PINKNOISE2(PROT, S)
%
%  Differs from PINKNOISE in that edges of spatial and 
%  temporal filters smoothed with Gaussian and Hamming 
%  windows.
%
% PROT is struct describing stimuli to be generated
% S is SCREEN structure
%
% see READPROT, SCREEN
% see PINKNOISE, PINKRND, PINKRND2

refreshRate = s.RefreshRate/s.DecimationRatio;

% each protocol consists of multiple stimuli
for iStim = 1:prot.nstim  
	sel = find([prot.pars.n]==iStim);
    
    % each stimulus consists of multiple segments
    for iSeg= 1:length(sel);
        fprintf(1,'Stimulus %i\n',iStim);
        pars = prot.pars(sel(iSeg));
        
        r = pinkrnd2(pars,refreshRate);

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


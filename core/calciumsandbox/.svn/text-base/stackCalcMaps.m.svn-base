function [maps]=stackCalcMaps(stack,stims,blanks, sigma, blankwin, stimwin)
%STACKCALCMAPS
% [MAPS]=STACKCALCMAPS(STACK,STIMS,BLANKS, SIGMA) where STACK is data, 
%  STIMS and BLANKS are cell arrays of size [NTRIALS,NSTIM] containing frame 
%  indices for the stimulation and blank epochs (see GETEPOCHS), 
%  SIGMA is width of smoothing Gaussian
% [MAPS]=STACKCALCMAPS(STACK,STIMS,BLANKS, SIGMA, BLANKWIN) specifies the
% number of frames to use during the blank epochs.
% [MAPS]=STACKCALCMAPS(STACK,STIMS,BLANKS, SIGMA, BLANKWIN, STIMWIN)
% specifies the numbers of frames during the stimulus epochs.
%
% set GETEPOCHS, MAPSCALCORI
%
% by Vincent Bonin

if nargin < 5
    blankwin = floor(getEpochsLen(blanks) / 2);
end

if nargin < 6 
    stimwin = getEpochsLen(stims);
end

g=fspecial('gaussian', [1 1]*sigma*5, sigma);    % spatial lowpass filter in pixel unit (after binning)
%g = 1;

% fprintf(1,'Calculating baseline...\n');
% 
% blankfr = [blanks{:,:}];
% maps.base = mean(stack(:,:,blankfr),3);

[ny,nx,nframes]=size(stack);

[ntrials,nstim]=size(blanks);

fprintf(1,'Calculating single trials...\n');

maps.trials.F = zeros(ny,nx,nstim,ntrials);
maps.trials.dF = zeros(ny,nx,nstim,ntrials);
maps.trials.ratio = zeros(ny,nx,nstim,ntrials);
maps.trials.basesd = zeros(ny,nx,nstim,ntrials);
maps.trials.zscore = zeros(ny,nx,nstim,ntrials);

for itrial = 1:ntrials
	fprintf(1,'Trial %i ',itrial);
    for istim = 1:nstim
        fprintf(1,'.');
        stimfr = [stims{itrial,istim}];stimfr = stimfr(end-stimwin+1:end);
        blankfr = [blanks{itrial,istim}];blankfr = blankfr(end-blankwin+1:end);
        
        maps.trials.F(:,:,istim,itrial) = filter2(g,mean(stack(:,:,stimfr),3));
        maps.trials.base(:,:,istim,itrial) = filter2(g,mean(stack(:,:,blankfr),3));
        maps.trials.basesd(:,:,istim,itrial) = filter2(g,std(double(stack(:,:,blankfr)),[],3));
        maps.trials.dF(:,:,istim,itrial) = maps.trials.F(:,:,istim,itrial)-maps.trials.base(:,:,istim,itrial);
        maps.trials.ratio(:,:,istim,itrial) = maps.trials.dF(:,:,istim,itrial)./maps.trials.base(:,:,istim,itrial)*100;
        maps.trials.zscore(:,:,istim,itrial) = maps.trials.dF(:,:,istim,itrial)./maps.trials.basesd(:,:,istim,itrial);
    end
    fprintf(1,'\n');    
end

fprintf(1,'Calculating trial averages...\n');

trials = 1:ntrials;

maps.av.F = zeros(ny,nx,nstim);
maps.av.dF = zeros(ny,nx,nstim);
maps.av.ratio = zeros(ny,nx,nstim);
maps.av.basesd = zeros(ny,nx,nstim);
maps.av.zscore = zeros(ny,nx,nstim);

for istim = 1:nstim
    fprintf(1,'Stimulus %i \n',istim);
    stimfr = [stims{trials,istim}];
    blankfr = [blanks{trials,istim}];
    maps.av.F(:,:,istim) = mean(maps.trials.F(:,:,istim,:),4);
    maps.av.base(:,:,istim) = mean(maps.trials.base(:,:,istim,:),4);    
    maps.av.dF(:,:,istim) = mean(maps.trials.dF(:,:,istim,:),4);
%    maps.av.ratio(:,:,istim) = mean(maps.trials.ratio(:,:,istim,:),4);
     maps.av.ratio(:,:,istim) = maps.av.dF(:,:,istim)./maps.av.F(:,:,istim)*100.0;
    maps.av.zscore(:,:,istim) = mean(maps.trials.zscore(:,:,istim,:),4);
end

fprintf(1,'Done!\n');

return;

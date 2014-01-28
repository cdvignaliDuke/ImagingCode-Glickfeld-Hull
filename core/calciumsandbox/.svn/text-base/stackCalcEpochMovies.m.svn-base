function [movies] = stackCalcEpochMovies(stack,epochs,blanks);
%STACKCALCEPOCHMOVIES Spits out bunch of movies as uint16
% [MOVIES] = STACKCALCEPOCHMOVIES(STACK,EPOCHS,BLANKS);
%
% code optimized for memory usage, modify with care
% still uses way too much memory. no reason why should not be output uint8

[ny,nx,nframes]=size(stack);

%% calculate bunch of movies.av
[nTrials,nStim]=size(epochs);

len = getEpochsLen(epochs);

movies.trials.dF = zeros(ny,nx,len,nTrials,nStim,'uint16');
movies.trials.ratio = zeros(ny,nx,len,nTrials,nStim,'uint16');

% single need for ratio calculation, a pain with integers
segment = zeros(ny,nx,len,'single');
baseline = zeros(ny,nx,'single');

offset = 2^15;

for iStim = 1:nStim
    fprintf(1,'Stimulus %i ',iStim);
    for iTrial = 1:nTrials
        fprintf(1,'.');

        ind = epochs{iTrial,iStim};      
        baseind = blanks{iTrial,iStim};
        
        ind = epochs{iTrial,iStim};
        
        baseind = blanks{iTrial,iStim};

        segment = single(stack(:,:,ind(1:len)));
        baseline = mean(stack(:,:,baseind),3);
        
        baseline(find(baseline(:)<1))=1;
        
        % calculates dF
        temp1 = bsxfun(@minus,segment,baseline);
        
        % calculates ratio
        temp2 = bsxfun(@rdivide,temp1,baseline)*100;
        
        % casts to uint16
        movies.trials.dF(:,:,:,iTrial,iStim) = uint16(temp1+offset);  
        movies.trials.ratio(:,:,:,iTrial,iStim) = uint16(temp2+offset);
    end    

    fprintf(1,'\n');
    
end

movies.av.dF = uint16(squeeze(mean(single(movies.trials.dF),4)));
movies.av.ratio = uint16(squeeze(mean(single(movies.trials.dF),4)));

return;

% function [movies] = stackCalcEpochMovies(stack,epochs,len, blankwin);
% % [movies] = stackCalcEpochMovies(stack,epochs blankwin);
% 
% if nargin < 3 
%     len = 200;
% end
% 
% if nargin < 4
%     blankwin = 64;
% end
% 
% [ny,nx,nframes]=size(stack);
% 
% %% calculate bunch of movies.av
% [nTrials,nStim]=size(epochs);
% 
% temp.trials.dF = zeros(ny,nx,len,nTrials,nStim,'single');
% temp.trials.ratio = zeros(ny,nx,len,nTrials,nStim,'single');
% 
% segment = zeros(ny,nx,len,'single');
% baseline = zeros(ny,nx,'single');
% 
% offset = int16(2^14);
% for iStim = 1:nStim
%     fprintf(1,'Stimulus %i ',iStim);
%     for iTrial = 1:nTrials
%         fprintf(1,'.');
%         ind = epochs{iTrial,iStim};
%         
%         if ~blankwin
%             baseind = ind(end-blankwin+1:end);
%         else
%             baseind = ind;
%         end
%         
%         segment = single(stack(:,:,ind(1:len)));
%         
%         baseline = mean(stack(:,:,baseind),3);
%         
%         temp.trials.dF(:,:,:,iTrial,iStim) = bsxfun(@minus,segment,baseline);
%         temp.trials.ratio(:,:,:,iTrial,iStim) = bsxfun(@rdivide,...
%                           temp.trials.dF(:,:,:,iTrial,iStim),baseline)*100;
%     end    
% 
%     fprintf(1,'\n');
%     
% end
% 
% temp.av.dF = squeeze(mean(temp.trials.dF,4));
% temp.av.ratio = squeeze(mean(temp.trials.dF,4));
% 
% movies.trials.dF = uint16(temp.trials.dF+2^15);
% movies.trials.ratio = uint16(temp.trials.ratio+2^15);
% 
% movies.av.dF = uint16(temp.av.dF+2^15);
% movies.av.ratio = uint16(temp.av.ratio+2^15);
% 
% return;
% 

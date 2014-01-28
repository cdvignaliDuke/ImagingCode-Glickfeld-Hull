function [stims,blanks,epochs]=getPresentationEpochs(p)
% [STIMS,BLANKS,EPOCHS]=GETPRESENTATIONEPOCHS(P)
% 
% see readPresentationLog, parsePresentationLog.

%%
nFramePulses  = length(p.FrameTimes);
nPictures = length(p.PictureTimes);

%%
nPicsPerTrial = p.info.nPicsPerStim * p.info.nStim;

% based on invalid assumption
%nExpectedPictures = p.info.nPicsPerStim * p.info.nStim * p.info.nTrials;

% fprintf(1,'Imaged Frames = %i. Frame Pulses = %i.\n',...
%            p.info.nImagedFrames,nFramePulses);
       
fprintf(1,'nPictures = %i.\n',nPictures);
%        
% if nExpectedPictures ~= nPictures
%     error('something is wrong');
% end

for iTrial = 1:p.info.nTrials
    for iStim = 1:p.info.nStim
        
        epochs{iTrial,iStim} = find(p.FrameStims==iStim & p.FrameTrials==iTrial);
        
        stims{iTrial,iStim} = find(p.FrameStims==iStim & p.FrameTrials==iTrial & p.FrameBlanks==0);
        
        blanks{iTrial,iStim} = find(p.FrameStims==iStim & p.FrameTrials==iTrial & p.FrameBlanks==1);        
        
    end
end

% lastframes = max(epochs{iTrial,iStim})+1:p.info.nImagedFrames;
% 
% switch p.info.blankMode
%     case 0
%         stims{iTrial,iStim} = [stims{iTrial,iStim} lastframes];
%         epochs{iTrial,iStim} = [epochs{iTrial,iStim} lastframes];
%     
%     case 1
%         stims{iTrial,iStim} = [stims{iTrial,iStim} lastframes];
%         epochs{iTrial,iStim} = [epochs{iTrial,iStim} lastframes];
%         
%     case 2        
%         blanks{iTrial,iStim} = [blanks{iTrial,iStim} lastframes];
%         epochs{iTrial,iStim} = [epochs{iTrial,iStim} lastframes];
%         
% end

frameSeqLen = length([epochs{:}]);

fprintf(1,'Frame Pulses = %i. Frame Seq Length = %i.\n',...
          nFramePulses,frameSeqLen);

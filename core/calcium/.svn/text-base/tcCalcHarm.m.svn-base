function harms = tcCalcHarm(timeCourses,Fs,F0)
%TCCALCHARM Computes harmonic responses
% HARMS = TCCALCHARM(TIMECOURSES,FS,F0) where TIMECOURSES is sorted time
% course array of size [NSAMPLES,NCELLS,NTRIALS,NSTIMS], FS is the sampling
% frequency, and F0 is the stimulus frequency.
%
% see TCSORTTRIALS
 
[nSamples,nCells,nTrials,nStims]= size(timeCourses);
tt = [0:nSamples-1]/Fs;
xx = exp(-j*2*pi*F0*tt(:))*ones(1,nStims);

harms =zeros(nCells,nTrials,nStims);

for iCell = 1: nCells
    for iTrial = 1:nTrials
        this = squeeze(timeCourses(:,iCell,iTrial,:));
        harms(iCell,iTrial,:)=sum(xx.*this)/Fs/tt(end);       
    end
end

if F0 > 0 
    harms = 4*harms;
end

return;

xx = [0:255]/frGetFrameRate;
test = sin(2*pi*xx)+2;
abs(tcCalcHarm(test(:),Fs,2))

%%
plot(tc)
tc = sum(timeCourses(:,1,1),3);
clf
xmax =max(tc);

vals =  .05:.05:4;
temp = [];
for f0 = vals
xx = exp(-j*2*pi*f0*tt(:))*ones(1,nStims);
temp(end+1)=abs(sum(tc.*xx));
end

plot(vals,temp)

%%
plot(tc.*xx,'k-o','markerfacecolor','w')
hold on;
plot(xx(1:32)*xmax,'k:')
xlim([-1 1]*xmax)
ylim([-1 1]*xmax);
grid on;
pause; 
end
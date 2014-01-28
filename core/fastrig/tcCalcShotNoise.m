function [snr,counts,slope,sd] = tcCalcShotNoise(timeCourses,w);
%TCCALCSHOTNOISE
%[SNR,COUNTS,SLOPE,SD] = TCCALCSHOTNOISE(TIMECOURSES,W);

if nargin < 2
    w = 512;
end

[nSamples,nCells]=size(timeCourses);

nSeg = floor(nSamples/w);
len = nSeg*w;

chopped = reshape(timeCourses(1:len,:),[w,nSeg,nCells]);

av = (squeeze(mean(chopped)));
mav = mean(av);

va = (squeeze(var(chopped)));

% shot noise imposes a lower bound on variance / mean ratio
slope = prctile(va(:)./av(:),5); % digital units / photon
% 
% clf;plot([0 1e3],[0 1e3],'k');hold on;plot(av,(va),'o')

% estimated shot noise variance
counts = mav ./ slope; % photons / bin
snr = sqrt(counts);
sd = sqrt(mav*slope);

return;

% up = 95;
% avmax = prctile(av,up);
% vamax = prctile(va,up);
% 
% avnorm = av(:)/slope;
% vanorm = va(:)/slope^2;
% 
% avnormmax = prctile(avnorm,up);
% vanormmax = prctile(vanorm,up);
% 
% % scatter plots
% figure;
% plot(avnorm(:),vanorm(:),'k.','markersize',5);
% hold on
% plot([0 avmax]/slope,f(slope,[0 avmax])/slope^2,'r','linewidth',2);
% xlim([0 avmax]/slope);ylim([0 vamax]/slope^2);
% 
% axis square;grid;
% xlabel('Mean (counts / sample)');
% ylabel('Variance');

return;
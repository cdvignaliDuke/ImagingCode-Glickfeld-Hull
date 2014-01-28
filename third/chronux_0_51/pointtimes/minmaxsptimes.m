function [mintime, maxtime]=minmaxsptimes(data)
% Find the minimum and maximum of the spike times in each channel
% Usage: [mintime, maxtime]=minmaxsptimes(data)
% Input:
% data  (spike times as a structural array of channels/trials dimensions; can also accept a 1d 
%               of spike times)
% Output:
% mintime       (minimum of the spike time across channels)
% maxtime       (maximum of the spike time across channels)
%
if isstruct(data)
   C=length(data);
   fnames=fieldnames(data);
   for ch=1:C
     eval(['dtmp=data(ch).' fnames{1} ';'])
     if ~isempty(dtmp)
        maxtime(ch)=max(dtmp);
        mintime(ch)=min(dtmp);
     else
        mintime(ch)=NaN;
        maxtime(ch)=NaN;
     end
   end;
   maxtime=max(maxtime); % maximum time
   mintime=min(mintime); % minimum time
else
     dtmp=data;
     if ~isempty(dtmp)
        maxtime=max(dtmp);
        mintime=min(dtmp);
     else
        mintime=NaN;
        maxtime=NaN;
     end
end
if mintime < 0 
   error('Minimum spike time is negative'); 
end
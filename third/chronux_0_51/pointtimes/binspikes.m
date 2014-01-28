function [dN,t]=binspikes(data,Fs,t)
% bin spikes at a specified frequency sampling i.e. sampling rate 1/sampling
% eg: 1ms accuracy use sampling = 1000
% Usage: [dN,t]=binspikes(data,Fs,t)
% Inputs:
% data   (data as a structure array of spike times; or as a single column
%        vector of spike times)
% Fs     (binning frequency)
% t      (the minimum and maximum times to be used to form the bins - [mint maxt]
%            - optional. Default use the spike times themselves to
%              determine the location of the bins. 
% Note: the times in data can be in any units. However, it is important
% that all units are chosen consistently. So, if spike times are in secs,
% Fs and t (if present) have to be in Hz and secs respectively. If spike
% times are in number of samples, Fs has to be 1, and t has to be in number
% of samples.
% Outputs:
% dN     (output binned spike counts as a matrix defined on bins starting with the
%         earliest spike across all channels and ending with the latest spike)
% t      (lower limit of each bin)
if nargin < 2; error('Need at least two input arguments'); end;
dt=1/Fs;
if isstruct(data);
   C=length(data);
   fnames=fieldnames(data);
   if nargin <3 || isempty(t);
       for ch=1:C
         eval(['dtmp=data(ch).' fnames{1} ';'])
         mintime(ch)=min(dtmp);
         maxtime(ch)=max(dtmp);
       end
       mintime=min(mintime);
       maxtime=max(maxtime);
       t=mintime:dt:maxtime;
   else
       for ch=1:C
         eval(['dtmp=data(ch).' fnames{1} ';'])
%          mintimech(ch)=min(dtmp);
         maxtimech(ch)=max(dtmp);
       end
       mintime=t(1);
       maxtime=t(end);
%        mintimech=min(mintimech);
       maxtimech=max(maxtimech);
       t=mintime:dt:maxtime;
       if maxtimech > max(t); t=[t maxtimech+dt]; end;
   end
   for ch=1:C;
       eval(['dtmp=data(ch).' fnames{1} ';'])
       x=histc(dtmp,t);
       dN(:,ch)=x(:);
   end
else
   dtmp=data;
   if nargin < 3;
      mintime=min(dtmp);
      maxtime=max(dtmp);
   else
      mintime=t(1);
      maxtime=t(end);
   end
   t=mintime:dt:maxtime;
   if max(dtmp)>max(t); t=[t maxtime+dt]; end;
   x=histc(dtmp,t);
   dN=x(:);
end
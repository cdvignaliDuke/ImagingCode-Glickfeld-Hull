function [raw]=formatframe(data,linestart,cyclelen,delay)

% maxdelay = 50;
% maxdelay = 100;
% 
% maxdelay=100; %SY

linelen = uint32(cyclelen/2);

ncycles = length(linestart);
% temp1 = chopvec(data,linestart+maxdelay-delay,linelen);
% temp2 = chopvec(data,linestart+linelen+maxdelay-delay,linelen);
temp1 = chopvec(data,linestart-delay,linelen);
temp2 = chopvec(data,linestart+linelen-delay,linelen);
temp3 = flipud(temp2);
temp4 = cat(1,temp1,temp3);
raw = reshape(temp4,[cyclelen/2,ncycles*2])';

% raw = cat(1,temp1,temp3);
% raw = reshape(raw,[cyclelen/2,ncycles*2])';


return;

% figure; imagesc(temp1); figure;imagesc(temp3);
function [delay,errs]=calculatedelay(data,linestart,cyclelen, delay, n)

persistent temp3;
persistent temp4; % temp1 and temp2 allocated in mex file

maxdelay = 50;
maxdelay=100; %SY
linelen = uint32(cyclelen/2);
ncycles = length(linestart);


if isnan(delay) && n > 0;
    delays =1:40;
else
    delays = delay-n:delay+n;
end

%%
errs = zeros(1,length(delays));
% tic;
% for iter = 1:100;
for idelay = 1:length(delays);
    temp1 = chopvec(data,linestart+maxdelay-delays(idelay),linelen);
    temp2 = chopvec(data,linestart+linelen+maxdelay-delays(idelay),linelen);
    temp3 = flipud(temp2);
    temp4 = abs(temp1-temp3); 
    errs(idelay) = sum(sum(temp4))/linelen/ncycles;
end
% end
% toc

[val,ind]=min(errs);
delay = delays(ind);
%err = errs(ind);


if ind==1 || ind==2 || ind==length(delays) || ind==length(delays)-1
    
[delay,errs]=calculatedelay(data,linestart,cyclelen, delay, n);   
    
end


%%
return;

% figure; imagesc(temp1); figure;imagesc(temp3);

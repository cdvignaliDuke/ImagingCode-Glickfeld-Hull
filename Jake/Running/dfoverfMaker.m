%% cell1 converts 60 frame averages into 60 frame df/f using forloop
f = mean(zreg(:,:,1:50),3);
n = 1;
dfoverf = zreg;
for n = 1:1500
    dfoverf(:,:,n) = (f - zreg(:,:,n))/f;
    n = n + size1;
end
clear f n
'done'

%% cell2  
%Alternate method using bsfxfun
f = mean(avgConValue1(:,:,1:7),3);
df = bsxfun(@minus,avgConValue1,f);
dfoverf = bsxfun(@rdivide,df,f);
clear f df
'done'
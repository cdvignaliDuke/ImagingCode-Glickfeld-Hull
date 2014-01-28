function b = tcGetBaseline(x)
%B = TCGETBASELINE(X)

fprintf('Calculating baseline ');

[nt,ncells]=size(x);

m = mean(x);
xm = bsxfun(@minus,x,m);

% cumulative sum to get F trajectory
xint = cumsum(xm);

ds = 16;
xint2 = tcDecimate(xint,ds);
xint3 = tcDecimate(xint,ds);

N = 4;                 % Order of polynomial fit
F = 21;                % Window length
[x0,x1,x2]=sgolaydiff(xint3,N,F);

b=[];
for icell = 1:ncells
    fprintf('.');
    % find breakpoints
    bp = crossing(x1(:,icell));
    
    % fit piecewise linear model
    xlin = x0(:,icell) - detrend(x0(:,icell),'linear',bp);
    
    xlinp = diff(xlin)/ds;
        
    ind = find(xlinp<0);
    
    %b(icell)=median(xlinp(ind));
    b(icell)=prctile(xlinp(ind),10);
    
%     clf;
%     plot(x(:,icell));
%     hold on;
%     plot(xlim,[1 1 ]*(b(icell)+m(icell)),'r')
%     pause(.5);

end
fprintf(' Done!\n');

b = b + m;

return;


%%
icell = 1;
x = (cumsum(tcremovedc(tcs2.correct(:,icell))));
k = gausswin(256);k = k/sum(k);
y = (filter(k,1,x));
yp = diff(y);
plot(yp);plot(xlinp(:,1)) 
bp = crossing(yp);
a = y - detrend(y,'linear',bp);
plot(diff(a))
% plot(y);
% hold on;
% plot(a,'r');

b = (diff(a).*(diff(a)<0));
b(find(~b))=nan;
c = nanmedian(b)
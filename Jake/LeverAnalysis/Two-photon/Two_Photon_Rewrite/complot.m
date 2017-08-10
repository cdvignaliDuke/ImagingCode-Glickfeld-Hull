function complot(sig, ICuse, dt)

for i = 1:length(ICuse)
    zsig(i, :) = zscore(sig(ICuse(i),:));
end

alpha = mean(max(zsig')-min(zsig'));
if islogical(zsig)
    alpha = 1.5*alpha;
end

zsig2 = zsig;
for i = 1:size(ICuse,2)
    zsig2(i,:) = zsig(i,:) - alpha*(i-1)*ones(size(zsig(1,:)));
end

tvec = [1:size(zsig,2)]*dt;
if islogical(zsig)
    plot(tvec, zsig2','LineWidth',1)
else
    plot(tvec, zsig2','LineWidth',1)
end
axis tight

set(gca,'YTick',(-size(zsig,1)+1)*alpha:alpha:0);
set(gca,'YTicklabel',fliplr(ICuse));

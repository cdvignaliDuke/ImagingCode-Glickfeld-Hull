
n = 10;
nz = size(img_reg,3);
img_me = zeros(npw,nph,nz/n);
for i = 1:nz/n
    i
    img_me(:,:,i) =  mean(img_reg(:,:, (i-1)*10+1:i*10),3);
end

cc = [];
for i = 1:5000
    f = floor((i-1)/10)+1;
    cc = [cc; corr2(img_reg(:,:,i), img_me(:,:,f))];
end

max_corr = squeeze(mean(mean(cc)));
mean(max_corr) - std(max_corr)/10

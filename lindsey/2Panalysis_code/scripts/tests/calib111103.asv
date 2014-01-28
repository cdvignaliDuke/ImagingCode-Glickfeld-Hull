fn = 'G:\users\lindsey\dataLG\LG111109_monitor\1on';
on1 = readtiff(fn);
siz = size(on1);

cyc = 20;

on1_tc = squeeze(mean(mean(on1(:,:,4:end),1),2));

ncyc = siz(3)/cyc;
on1_rep_tc= zeros(cyc,ncyc);
start = 5;
for icyc = 1:ncyc-1
    on1_rep_tc(:,icyc) = on1_tc(start:start+cyc-1,:);
    start = start+cyc;
end

on1_rep_tc_dF = bsxfun(@minus, on1_rep_tc, mean(on1_rep_tc,1));
figure;
plot(on1_rep_tc_dF); 
figure;
plot(on16_rep_tc(:,1:83))
hold on 
plot(on16_rep_tc(:,85))

[min_tc ind] = min(min(on2_rep_tc,[],2),[],1);
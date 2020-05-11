function [reg shift] = shift_opt(data, target, n);
% finds optimal shift to align data and target using correlation 
% coefficient with pixel steps n in each direction
orig_r = triu2vec(corrcoef(data, target));
sz = size(data);
range = [-n:n];
out_n = length(range).^2;
out_grid = zeros(out_n,4);
reg_temp = zeros(sz(1),sz(2),out_n);
r = zeros(1,out_n);
start = 1;
for i = 1:length(range);
    out_grid(start:start+length(range)-1,3) = range(i);
    for j = 1:length(range);
        out_grid(start,4) = range(j);
        [o reg_temp(:,:,start)] = stackRegister_MA(data, [], [], out_grid(start,:));
        r(start) = triu2vec(corrcoef(reg_temp(:,:,start), target));
        start = start+1;
    end
end
        
[r_max r_ind] = max(r,[],2);
if r_max>orig_r
    reg = reg_temp(:,:,r_ind);
    shift = out_grid(r_ind,:);
else
    reg = data;
    shift = out_grid(61,:);
end
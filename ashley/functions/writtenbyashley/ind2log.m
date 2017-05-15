function idx_log = ind2log(ind_sub,n)
    
    x = zeros(1,n);
    x(ind_sub) = 1;
    x = logical(x);
    
    idx_log = x;
end
function err = ste(array, dim)

    stdev = nanstd(array,[],dim);
    n = size(array,dim);
    
    err = stdev./(sqrt(n-1));

end
function err = ste(array, dim)

    stdev = std(array,[],dim);
    n = size(array,dim);
    
    err = stdev./(sqrt(n-1));

end
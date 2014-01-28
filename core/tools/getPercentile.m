function value=getPercentile(any_matrix, percent)

a=sort(any_matrix(:));
value=a(round(length(a)*percent/100));
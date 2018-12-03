function si = getSelectivityIndex(responses1,responses2)

avg1 = mean(responses1,1);
avg2 = mean(responses2,1);
var1 = var(responses1,[],1);
var2 = var(responses2,[],1);
n1 = size(responses1,1)-1;
n2 = size(responses2,1)-1;

stdpool = sqrt((var1*n1) + (var2*n2)) ./ (n1 + n2);

si = (avg1 - avg2) ./ stdpool;

end

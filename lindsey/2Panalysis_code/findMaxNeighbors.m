function out = findMaxNeighbors(X,D)
%finds point in X that is within D (distance) from most other X
for i = 1:size(X,2)
    dists = sqrt((X(1,:)-X(1,i)).^2 + (X(2,:)-X(2,i)).^2);
    n = length(find(dists<D));
    if i == 1
        n_max = n;
        out = X(:,i);
    elseif n > n_max
        n_max = n;
        out = X(:,i);
    end
end
        
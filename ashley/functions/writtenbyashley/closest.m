function [minVal Y] = closest(X,VAL)
% find the index of the value in a vector, X, that is closest to indicated value, VAL
tmp = abs(X-VAL);
[minAbs Y]  = min(tmp);
minVal = X(Y);

function y = fano(x,dim)
%fano factor calculation, variance/mean
%   For vectors, Y = FANO(X) returns the fano factor.  For matrices,
%   Y is a row vector containing the fano factor of each column.  
%   
%   Y = FANO(X,DIM) returns the fano factor of X along dimension DIM.
if nargin < 2
    y = std(x)./mean(x);
else
    y = std(x,[],dim)./mean(x,dim);
end
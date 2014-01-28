function y = VonMises2(A, xdata)

% double von Mises function to fit direction tuning curve
% A: paramters of von Mises function
%       A(1) = A1, A(2)=A2, A(3) = k1, A(4) = k1, A(5) = phi1, A(6) = phi2
%       y = A1 * exp( k1 * cos (x-phi1) - 1) + A2 * exp(k2 * cos (x-phi2) -1) 
% xdata: a vector which contains x values in degree unit (not radian). 
%
%   Kenichi Ohki 09/19/04
%

y = A(1) * exp( A(3) * (cos((xdata-A(5))*pi/180) -1 ))...
    + A(2) * exp( A(4) * (cos((xdata-A(6))*pi/180) -1 )) ;

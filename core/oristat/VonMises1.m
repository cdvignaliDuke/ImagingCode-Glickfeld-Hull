function y = VonMises1(A, xdata)

% double von Mises function to fit direction tuning curve
% A: paramters of von Mises function
%       A(1) = A, A(2) = k, A(3) = phi
%       y = A * exp( k * cos 2 * (x-phi) - 1) 
% xdata: a vector which contains x values in degree unit (not radian). 
%
%   Kenichi Ohki 09/19/04
%

y = A(1) * exp( A(2) * (cos((xdata-A(3))*pi/90) -1 ));

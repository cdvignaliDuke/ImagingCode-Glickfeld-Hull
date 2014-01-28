function [yi] = interp0SB(x,y,xi)
% interp0SB: zero order interpolation function (lookup table)
% If of limits then the extreme points in y are taken as output. 
% interp0SB can be used together with MEX simulation files. For MEX
% simulation functions it is IMPORTANT that the elements of the x and y
% vectors are numeric and SEPARATED BY COMMATA! 
% 
% USAGE:
% ======
% [yi] = interp0SB(x,y,xi)   
%
% x: vector of function arguments
% y: vector of function values at the points given by x
% xi: scalar value for which to determine y by zero order interpolation

% Information:
% ============
% Copyright (C) 2008 Henning Schmidt, henning@sbtoolbox2.org
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
% USA.

yi = zeros(1,length(xi));
for k=1:length(xi),
    if xi(k) < x(2),
        yi(k) = y(1);
    elseif xi(k) >= x(end),
        yi(k) = y(end);
    else
        yi(k) = y(max(find(x<=xi(k))));
    end
end
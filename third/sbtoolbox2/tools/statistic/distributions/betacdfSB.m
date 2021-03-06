function cdf = betacdfSB(x, a, b)
% Cumulative density function of the Beta distribution
%
% USAGE:
% ======
% cdf = betacdfSB(x, a, b)
%
% For each element of 'x', returns the CDF at 'x' of the beta
% distribution with parameters 'a' and 'b', i.e.,
% PROB (beta ('a', 'b') <= 'x').

% Information:
% ============
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>
% This function has been taken from Octave and adapted for the SBTOOLBOX2
% by Henning Schmidt
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

  if (nargin ~= 3)
    error('Incorrect number of input arguments'); 
  end

  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('x, a and b must be of common size or scalar');
    end
  end

  sz = size(x);
  cdf = zeros (sz);

  k = find (~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    cdf (k) = NaN;
  end

  k = find ((x >= 1) & (a > 0) & (b > 0));
  if (any (k))
    cdf (k) = 1;
  end

  k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
  if (any (k))
    if (isscalar (a) && isscalar(b))
      cdf (k) = betainc(x(k), a, b);
    else
      cdf (k) = betainc(x(k), a(k), b(k));
    end
  end

end

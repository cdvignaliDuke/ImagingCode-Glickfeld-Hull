function [output] = quadprogSB(D, d, A, varargin)
% quadprogSB: Solve a Quadratic Programming Problem using the 
% dual method of Goldfarb and Idnani (1982, 1983).
% 
%     min(-d'*x + 1/2*x'*D*x)  
%      x
%
%     subject to 
%                  A x <= b
%
% USAGE:
% ======
% [output] = quadprogSB(D, d, A)
% [output] = quadprogSB(D, d, A, b)
% [output] = quadprogSB(D, d, A, b, meq)
% [output] = quadprogSB(D, d, A, b, meq, factorized)
%
% D:     matrix appearing in the quadratic function to be minimized.  
% d:     vector appearing in the quadratic function to be minimized.  
% A:     matrix defining the constraints under which we want to minimize
%        the quadratic function. 
% b:     vector holding the values of b_0 (defaults to zero).  
% meq:   the first meq constraints are treated as equality constraints, 
%        all further as inequality constraints (defaults to 0).  
% factorized: if 1, then we are passing R^(-1) (where D = R^T R) 
%        instead of the matrix D in the argument Dmat. 
%
% DEFAULT VALUES:
% ===============
% b:          zero vector
% meq:        0
% factorized: 0 
%
% Output Arguments:
% =================
% The output argument is a structure with the following field:
% 
% output.sol:      nx1 the final solution (x in the notation above)
% output.soluc:    contains on exit the solution to the initial, i.e.,
%                  unconstrained problem
% output.optval:   scalar, the value of the criterion at the minimum
% output.actconst: indices of the constraints which are active in the final
%                  fit (int)
% output.iter:     the number of iterations
% output.delconst: how many constraints were deleted after they became
%                  active 
%
% References:
% ===========
% Goldfarb, D. and Idnani, A. (1982). Dual and Primal-Dual Methods for
% Solving Strictly Convex Quadratic Programs. In Numerical Analysis J.P.
% Hennart, ed. Springer-Verlag, Berlin. pp. 226-239.   
%
% Goldfarb, D. and Idnani, A. (1983). A numerically stable dual method for
% solving strictly convex quadratic programs. Mathematical Programming 27,
% 1-33.   

% Information:
% ============
% Copyright (C) 2008  Henning Schmidt, henning@sbtoolbox2.org
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

% Check number of inputs
if nargin < 3 || nargin > 6, error('Incorrect number of input arguments.'); end

% Convert the inputs to the format required by the Fortran code
% (min(-d^T x + 1/2 x^T D x) with the constraints A^T x >= b_0)
Dmat = D;
dvec = -d(:);
Amat = -A';
% handle bvec further below

% Get dimensions
n = size(Dmat,1);
q = size(Amat,2);

% Assign defaults
bvec = zeros(1,q);
meq = 0;
factorized = 0;

% Handle variable input arguments
if nargin >= 4, bvec = varargin{1}; bvec = -bvec(:); end
if nargin >= 5, meq = varargin{2}; end
if nargin == 6, factorized = varargin{3}; end

% Error checks
if n ~= size(Dmat,2), error('Dmat is not symmetric!'); end
if n ~= length(dvec), error('Dmat and dvec are incompatible!'); end
if n ~= size(Amat,1), error('Amat and dvec are incompatible!'); end
if q ~= length(bvec), error('Amat and bvec are incompatible!'); end
if meq > q || meq < 0, error('Value of meq is invalid!'); end
  
% Calling Fortran second time
Dmatdummy = Dmat+1;         % This strange construct forces MATLAB to really
Dmatdummy = Dmatdummy-1;    % make Dmat local. Otherwise it changes the content
                            % even of the D variable used to call this
                            % function.
[sol,soluc,crval,iact,nact,iter,ierr] = solveQP(Dmatdummy,dvec,n,n,Amat,bvec,n,q,meq,factorized);

% Check result
if ierr == 1, error('Constraints are inconsistent, no solution!'); end
if ierr == 2, error('Matrix D in quadratic function is not positive definite!'); end

if nact == 0,
    actconst = [];
else
    actconst = iact(1:nact);
end

% Output handling
output = [];
output.sol = sol(:);
output.soluc = soluc(:);
output.optval = crval;
output.actconst = actconst;
output.iter = iter(1);
output.delconst = iter(2);

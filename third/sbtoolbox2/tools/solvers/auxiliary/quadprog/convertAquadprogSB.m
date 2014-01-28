function [Amat,Aind] = convertAquadprogSB(A)
% convertAquadprogSB: Converts a matrix A into the matrics Amat and Aind
% that are required by the quadprogcompactSB solver.
% 
% The matrix A is the one that is used to define the constraints for the
% quadratic programming solver quadprogSB.
%
%           A x <= b
%
% In case of a huge sparsely populated A matrix it is advantageous to use
% the compact solver quadprogcompactSB instead. The matrix A is then
% defined by the two matrices Amat and Aind that have the following
% meaning:
%
% Amat:  matrix containing the non-zero elements of the matrix A that 
%        defines the constraints. If m_i denotes the number of non-zero 
%        elements in the i-th column of THE TRANSPOSED matrix A then the
%        first m_i entries of the i-th column of Amat hold these non-zero
%        elements (AGAIN OF THE TRANSPOSED MATRIX A). (If maxmi denotes the
%        maximum of all m_i, then each column of Amat may have arbitrary
%        elements from row m_i+1 to row maxmi in the i-th column.)      
% Aind:  matrix of integers. The first element of each column gives the
%        number of non-zero elements in the corresponding column of the
%        TRANSPOSED matrix A. The following entries in each column contain
%        the indexes of the rows in which these non-zero elements are.     
%
% This function here can be used to convert A to Aind and Amat.

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

% Transpose A
Atrans = A';

% Get dimensions for Aind and Amat
ncola = size(Atrans,2);
nrowa = max(sum(Atrans~=0));

% Initialize Aind and Amat
Aind = zeros(nrowa+1,ncola);
Amat = zeros(nrowa,ncola);

Aind(1,:) = sum(Atrans~=0);
for k=1:size(Atrans,2),
    nonzeroind = find(Atrans(:,k)~=0);
    nonzero = Atrans(nonzeroind,k);
    Aind(2:length(nonzeroind)+1,k) = nonzeroind;
    Amat(1:length(nonzero),k) = nonzero;
end

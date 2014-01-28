%CPD_TRANSFORM Applies the transformation to point set.
%   X = CPD_TRANSFORM(Y, normal);
%   applies the transformation found by CPD to the Y data set. "normal" is
%   a structure of transformation parameters.
%
%   Input
%   ------------------ 
%   Y                  real NxD, full 2-D matrix of the point set to transform. 
%                       N - number of points, D- dimensions. 
%   normal          transformation parameters structure, 
%                       second output of CPD_REGISTER function
%
%   Output
%   ------------------ 
%   X              Transformed version of Y
%
%   Examples
%   --------
%    See cpd_example2.m
%
%   See also CPD_REGISTER.

% Copyright (C) 2006 Andriy Myronenko (myron@csee.ogi.edu)
%
%     This file is part of the Coherent Point Drift (CPD) package.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     CPD package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with CPD package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

function X = cpd_transform(Y, normal)

if nargin<2, error('cpd_transform error! Not enough input parameters.'); end;

% prescale to normal system
[m, d]=size(Y);
Y=Y;
Y=Y;

% Apply transformation
G=cpd_G(Y, normal.Y0, normal.beta);
X=Y+G*normal.W;

% scale back

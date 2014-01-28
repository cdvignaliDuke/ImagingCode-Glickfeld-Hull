%CPD_P Construct matrix P of the posterior probabilities
%   [P, E]=CPD_P(x, y, sigma, [outliers]) returns matrix P of the posterior
%   probabilities for the GMM and negative log-likelihood value E. if outliers >0 also assumes presence of
%   extra uniform distribution to account for the outliers.
%   The function returnes empty P, E if sigma value is too small, this means
%   it's time to stop the annealing procedure.
%
%   Input
%   ------------------ 
%   x          MxD real, full 2-D matrices of GGM centroids locations
%   y          NxD real, full 2-D matrices of data point locations
%   sigma      std of the gaussians.
%   outliers   inverse of the domain size of the uniform distribution. 0 - don't use any extra uniform distribution. 
%
%   Output
%   ------------------ 
%   P          Matrix of the posterior probabilities (MxN)
%   E          Negative log-likelihood value
%
%   Examples
%   --------
%       x= [1 2; 3 4; 5 6;];
%       y=x+1;
%       sigma=2;
%       outliers=0;
%       [P, E]=cpd_P(x, y, sigma, outliers);
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

function [P, E]=cpd_P(x, y, sigma, outliers)

if nargin<3, error('cpd_P.m error! Not enough input parameters.'); end;
if ~exist('outliers','var') || isempty(outliers), outliers = 0; end;

k=-2*sigma^2;
[n, d]=size(x);[m, d]=size(y);

P=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);

P=squeeze(sum(P.^2,2));
P=P/k;
P=exp(P);

if outliers
  Pn=outliers*(-k*pi)^(0.5*d)*ones(1,m);
  s=sum([P;Pn]);
else
  s=sum(P);
end
  
if nnz(s)==numel(s)
    E=-sum(log(s));
    P=P./repmat(s,n,1);
else
    P=[];E=[];
end
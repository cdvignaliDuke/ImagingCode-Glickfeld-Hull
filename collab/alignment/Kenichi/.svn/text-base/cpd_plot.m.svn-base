%CPD_PLOT Plot the CPD current iteration result
%   CPD_PLOT(X, Y0, Y, W, sigma, [disk, sampling]); plots the result on
%   current iteration. Works only for 2D and 3D data sets.
%
%   Input
%   ------------------ 
%   X           Reference point set matrix NxD;
%   Y0          Template point set matrix MxD (initial postions of GMM centroids);
%   Y           Current postions of GMM centroids;
%   W           Transformation parameters;
%   sigma       Std of the gaussians in the GMM
%   disk        Discretisation of correspondances lines. disk=1 show all lines.
%               disk=2 - every second, etc.
%   sampling    Sampling of the mesh to show the transformation.
%
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

function cpd_plot(X, Y0, Y, W, sigma, disk, sampling)

if nargin<5, error('cpd_plot.m error! Not enough input parameters.'); end;
if ~exist('disk','var') || isempty(disk), disk = 1; end;
if ~exist('sampling','var') || isempty(sampling), sampling = 0.1; end;

[m, d]=size(Y0);

if d>3, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;
if d<2, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;

% for 2D case
if d==2,

    subplot(2,2,1); plot(X(:,1), X(:,2),'r*', Y(:,1), Y(:,2),'bo'); title('X data (red). Y GMM centroids (blue)');
    subplot(2,2,2); plot( X(:,1), X(:,2),'r*', Y0(:,1), Y0(:,2),'bo' ); title('X data (red). Y0 initial positions of GMM centroids (blue)');
    subplot(2,2,3); plot(Y(:,1), Y(:,2),'bo', Y0(:,1), Y0(:,2),'mo' ); title('Dispacements of Y (GMM cetroids)');
    hold on;
    for i=1:disk:m, plot([Y(i,1) Y0(i,1)],[Y(i,2) Y0(i,2)],':'); end;
    hold off;

    %%% show the transformation field
    subplot(2,2,4);
    xymin=min(Y);
    xymax=max(Y);

    [xn, yn]=meshgrid(xymin(1):sampling:sampling+xymax(1),xymin(2):sampling:sampling+xymax(2));
    xy=[xn(:) yn(:)];
    G=cpd_G(xy,Y0, sigma);
    Y2=xy+G*W;

    [m, n] = size(xn);
    A = spdiags(ones(m*n,2),[1 m],m*n,m*n);
    a=m:m:size(xy,1)-1;
    A(sub2ind([m*n m*n],a,a+1))=0;

    gplot(A,xy,':');
    hold on;gplot(A,Y2);
    hold off;
    title('Deformation field from Y0 to Y');
   

else
% for 3D case
    subplot(1,2,1);  plot3(X(:,1),X(:,2),X(:,3),'r.');hold on;plot3(Y(:,1),Y(:,2),Y(:,3),'bo'); hold off;  title('X data (red). Y GMM centroids (blue)');set(gca,'CameraPosition',[15 -50 8]);
    subplot(1,2,2);  plot3(X(:,1),X(:,2),X(:,3),'r.');hold on;plot3(Y0(:,1),Y0(:,2),Y0(:,3),'bo'); hold off;  title('X data (red). Y0 initial positions of GMM centroids (blue)');set(gca,'CameraPosition',[15 -50 8]);
end

drawnow;
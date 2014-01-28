%% Apply affine transformation on findcell3D results (inputs: coordf, coordv, centroidsf, centroidsv; output: tformed, tformed 2)
% construct the affine matirx, M, from hand-picked coordinates
% input: hand-picked coordinates from the fixed and vivo image (coordf & coordv)
%        raw data of hand-picked coordinates from the fixed and vivo image (excell)
% output: 3-by-4 affine matrix, M 
load 'D:\ReidLab rotation\DAPI_stacks\test\excell';
load 'D:\ReidLab rotation\DAPI_stacks\test\BW_fixed_z134_z301\labeled3D_centroids.mat';
load 'D:\ReidLab rotation\DAPI_stacks\test\BW_Set2_z80_z367\labeled3D_centroids.mat';
excell(:,3)=excell(:,3)+1; % excell file (from imageJ) starts from z=0, while findcell3D (matlab code) starts from z=1
i=[1:2:size(excell,1)];
j=[2:2:size(excell,1)];
coordf=excell(i,:);coordf=(coordf)';
coordf=cat(1, coordf,ones(1,size(excell,1)/2));
coordv=excell(j,:); coordv=(coordv)';
M=coordv/coordf;
em=(coordv-M*coordf).^2; 
errors=sqrt(em(1,:)+em(2,:)+em(3,:)); % errors in pixel of transformed coordinates. 

tformed=M*centroidsf;
tformed2=M*coordf
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% Apply CPD transformation on affined resultes (input: centroidsv, tformed , coordv, tformed2; output: tformedCPDed, tformed2CPDed)

% set all coordinates in N*3 structure.
coordf=coordf'; coordv=coordv';
centroidsv=centroidsv'; centroidsf=centroidsf';
tformed2=tformed2'; tformed=tformed';

%  CPD optional parameters are all shown with default values.
%  You can also use just  Y=cpd_register(X, Y0);

viz=1;         % visualize each iteration
outliers=0;    % don't assume any outliers
sigma=1;       % start annealing from sigma=1;
beta=1;        % std of G. Smaller value allows local deformations. Large - almost rigid
lambda=1;      % weight of regularization
anneal=0.97;   % annealing constant
max_it=100;    % maximum number of iterations
tol=1e-5;      % stopping tolerance 
EMtol=1e-2;    % EM stopping tolerance

% CPD registration

%--------------------------------------------------------------------------
%------------------------------------------------------------------------

%CPD_REGISTER Non-rigid registration of two point sets.
%   [Y, normal]=CPD_REGISTER(X,Y0, [outliers, sigma, beta, lambda, anneal, max_it, tol, EMtol]);
%   returns the registered (aligned) version of Y (template point
%   set) over the X (reference point set) and parameters of the
%   transformation W.
%
%   Input
%   ------------------ 
%   X       real NxD, full 2-D matrix of reference point set. 
%            N - number of points, D- dimensions.  
%   Y0      real MxD, full 2-D matrix of template point set. 
%            M - number of points, D- dimensions.  
%
%  viz           viz=1 to visualize each iteration. (default 0) 
%  outliers      Outliers' strength. outliers = 0 for no outliers.
%                   outliers > 0 to represent prior knowledge on amount of outliers 
%  sigma        Initial value of std in GMM.
%  beta          Std of affinity matrix G. Controls the smoothness of the transformation
%  lambda      Strength of the regularization
%  anneal       Annealing rate [0.93 0.98]
%  max_it       Maximum number of iterations
%  tol             Tolerance
%  EMtol         EM tolerance
%
%   Output
%   ------------------ 
%   Y              Final postion of GMM centroids. Aligned final positions
%                   of template point set onto reference point set.
%   normal      Structer of transformation parameters that can be used to
%                   apply found transformation to any other data set
%
%   Examples
%   --------
%     X=[0 1; 1 0; 0 0; 1 1];
%     Y0=X+0.01*randn(size(X))+0.2;
%     Y = cpd_register(X, Y0);
% 
%     viz=1;
%     outliers=0;
%     sigma=1;
%     beta=1;
%     lambda=1;
%     anneal=0.97;
%     max_it=100;
%     tol=1e-5;
%     EMtol=1e-2;
%     [Y, normal] = cpd_register(X, Y0, viz, outliers, sigma, beta, lambda,anneal, max_it, tol, EMtol);


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



if nargin<2, error('CPD error! Not enough input parameters.'); end;
if ~exist('viz','var') || isempty(viz), viz=0; end;
if ~exist('outliers','var') || isempty(outliers), outliers = 0; end;
if ~exist('sigma','var') || isempty(sigma), sigma = 1; end;
if ~exist('beta','var') || isempty(beta), beta = 1; end;
if ~exist('lambda','var') || isempty(lambda), lambda = 1; end;
if ~exist('anneal','var') || isempty(anneal), anneal = 0.97; end;
if ~exist('max_it','var') || isempty(max_it), max_it = 150; end;
if ~exist('tol','var') || isempty(tol), tol = 1e-5; end;
if ~exist('EMtol','var') || isempty(EMtol), EMtol = 1e-3; end;


% Preprocessing and normalization

Y0orig=tformed2;X=double(coordv);Y0=double(tformed2);
%X=unique(X,'rows'); Y0=unique(Y0,'rows'); 
[n, d]=size(X); [m, d]=size(Y0);

% Rescaling and shifting to the origin
[X, Y0, normal]= cpd_normalize(X, Y0);
Y=Y0; iter=0; E=1; ntol=tol+10; W=zeros(m,d);

% Construct affinity matrix G
G=cpd_G(Y0,Y0,beta);

% Annealing. Itarate for different sigma.
while (iter<max_it) && (ntol > tol)
    EMiter=0; EMtol=tol+10; Yold=Y;
    
    % EM iterations for constant sigma.
    while (EMiter<max_it)  && (EMtol > tol)
        
        % E-step. Constuct posterior probability matrix P
        [P, Eu]=cpd_P(Y, X, sigma, outliers);
        
            % STOP if simga is already annealed to very small value.
            if isempty(P), iter=max_it; disp('Annealing process has finished'); break; end;

            E_old=E; E=Eu+lambda/2*trace(W'*G*W); % CPD energy function.
            disp(['CPD iter = ' num2str(iter) '  EMiter= ' num2str(EMiter) ' E= ' num2str(E) ' sigma= ' num2str(sigma)]);
                
        % M-step. Solve linear system for W.
        dP=spdiags(sum(P,2),0,m,m); % precompute diag(P)
        W=(dP*G+lambda*sigma^2*eye(m))\(P*X-dP*Y0);
        
            % update Y postions
            Y=Y0+G*W;
            EMtol=norm((E_old-E)/E_old);
            EMiter=EMiter+1;   
    end
    
    
    % Anneal
    sigma=sigma*anneal;                              
    iter=iter+1;
    ntol=norm(Yold-Y);
    
end

normal.W=W;
normal.Y0=Y0;
normal.beta=beta;
tformed2CPDed = cpd_transform(Y0orig, normal);

disp('CPD registration succesfully completed.');


%--------------------------------------------------------------------------
%-------------------------------------------------------------------------
% CPD registration for centroids
figure;
tformedCPDed=cpd_register(centroidsv, tformed, viz, outliers, sigma, beta, lambda, anneal, max_it, tol, EMtol);



% plot CPDed results. 
figure;scatter3(centroidsv(:,1),centroidsv(:,2),centroidsv(:,3),'o','r');
hold on; scatter3(coordv(:,1),coordv(:,2),coordv(:,3),'filled','g');
axis([0 600 0 600 -100 600]);
h=gcf;
set(h,'name','vivo2');

figure;scatter3(tformedCPDed(:,1),tformedCPDed(:,2),tformedCPDed(:,3),'o','b');
hold on;scatter3(tformed2CPDed(:,1),tformed2CPDed(:,2),tformed2CPDed(:,3),'filled','g');
axis([0 600 0 600 -100 600]);
h=gcf;
set(h,'name','fixed_tformed_CPDed');

% set all coordinates back to 3*N structure
centroidsf=centroidsf';
centroidsv=centroidsv';
coordf=coordf';
coordv=coordv';
tformed=tformed';
tformed2=tformed2';
tformed2CPDed=tformed2CPDed';
tformedCPDed=tformedCPDed';
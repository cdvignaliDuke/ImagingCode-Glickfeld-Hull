%VB Gabor fitting, modified to fit 2D elliptical 

%params VB sends in: sta.kmean

s = struct;
s.data = data;

%m = @(pars,dummy) Gauss2D_ellipseMA(pars,grid2.sfsf,grid2.tftf);
m = @(pars,sftf) Gauss2D_ellipseMA(pars,sftf);
s.m = func2str(m);

% A = pars(1);
% sigma_SF = pars(2); % in cycles/deg
% sigma_TF = pars(3); % in Hz
% sf0 = pars(4); %center sf
% tf0 = pars(5); %center tf

%s.lb =  [.001 1/6   .005  0 xmin xmin -pi/2];
%s.ub =  [1000 2/6   .05   2*pi  xmax  xmax  pi/2];

clear('x');
clear('y');
[m2,n2] = size(grid2.sfsf);
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));
y = s.data(:);

Mask_NaN = isnan(y);
ind_NaN = find(isnan(y));
ind_noNaN = find(Mask_NaN == 0);
if ~isempty(ind_NaN)
    x(ind_NaN,:) = [];
    y(ind_NaN) = [];
end


%A sigma_SF sigma_TF sf0 tf0 xi
%s.lb =  [.001 .25 .25  (min(x(:,1))) (min(x(:,2))) -2];
%s.ub =  [3 3   3    (max(x(:,1))) (max(x(:,2))) 2];

s.lb =  [.001 .25 .25  (min(x(:,1))) (min(x(:,2))) -2];
s.ub =  [3 4   4    (max(x(:,1))) (max(x(:,2))) 2];

%OLD WAS up to 110514:
%s.lb =  [.001 .25 .25  (min(x(:,1))-1) (min(x(:,2))-1) -2];
%s.ub =  [1000 3   3    (max(x(:,1))+1) (max(x(:,2))+1) 2];

%s.opt = 'silent';

%ph = 0:pi/2:3*pi/2;
%ori = 0:pi/4:3*pi/4;
%pos = [-15,-15;-15 0;-15 15;0 -15; 0 0;0 15; 15 -15; 15 0 ; 15 15];

Nsamps = 2;
dbin = [s.ub(1) - s.lb(1); ...
    s.ub(2) - s.lb(2); ...
    s.ub(3) - s.lb(3); ...
    s.ub(4) - s.lb(4); ...
    s.ub(5) - s.lb(5); ...
    s.ub(6) - s.lb(6)] ./ (Nsamps-1);
    
sigma_SF_vec = s.lb(2)+dbin(2)/2:dbin(2):s.ub(2);
sigma_TF_vec = s.lb(3)+dbin(3)/2:dbin(3):s.ub(3);
SF_vec = s.lb(4)+dbin(4)/2:dbin(4):s.ub(4);
TF_vec = s.lb(5)+dbin(5)/2:dbin(5):s.ub(5);
xi_vec = s.lb(6)+dbin(6)/2:dbin(6):s.ub(6);


clear temp;
index = 1;

% loop over initial parameters not to get stuck in local minimum

names = {'x2','resnorm'};
args = cell(5);

%options.Display = 'off';

for iSigSF = 1:length(sigma_SF_vec)
    for iSigTF = 1:length(sigma_TF_vec)
        for iSF = 1:length(SF_vec)
            for iTF = 1:length(TF_vec)
                for ixi = 1:length(xi_vec)
%                    fprintf('.');
                    s.x0 =  [max(max(s.data)) sigma_SF_vec(iSigSF) sigma_TF_vec(iSigTF) SF_vec(iSF) TF_vec(iTF) xi_vec(ixi)];

                    %                test = Gauss2D_ellipseMA(s.x0,x);

                    %                [args{:}] = lsqcurvefit(m,s.x0,grid2.sfsf,grid2.tftf,s.data,s.lb,s.ub,s.opt);
                    %                [args{:}] =
                    %                lsqcurvefit(m,s.x0,grid2.sfsf,grid2.tftf,s.data,s.lb,s.ub);
                    %                [args{:}] = lsqcurvefit(m,s.x0,x,y,s.lb,s.ub);
                    options = optimset('Display', 'off');
                    [x2,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]  = lsqcurvefit(m,s.x0,x,y,s.lb,s.ub,options);
                    temp(index) = cell2struct({x2,Resnorm},names,2);
                    index = index +1;
                end
            end;
        end;
    end;
end


[val,ind]=min([temp.resnorm]);
tmp = [];
for countres = 1:length(temp)
    tmp = [tmp; temp(countres).resnorm temp(countres).x2];
end

s.fit = temp(ind);
s.x = s.fit.x2;
s.k2 = zeros(m2*n2,1);
s.k2(ind_noNaN) = m(s.fit.x2,x);
s.k2b = reshape(s.k2,m2,n2);

s.res = s.data-s.k2b;
s.Maxfit = max(max(s.k2b));
s.Maxdata = max(max(s.data));


%find highcut for SF and TF 
%x(:,1) = log2(grid2.sfsf(:));
%x(:,2) = log2(grid2.tftf(:));
%    TF_vec0 = ind_TFuse(:,2);
%    SF_vec0 = ind_SFuse(:,2);
x00 = zeros(size(grid2.sfsf00,1)*size(grid2.sfsf00,2),2);
x00(:,1) = (grid2.sfsf00(:));
x00(:,2) = (grid2.tftf00(:));
k2b00 = m(s.fit.x2,x00);
%now find 10% and 50%
s.Maxfit00 = max(max(k2b00));

%tmp2 = max(tmp,[],1); %max over TFs
%tmp3 = tmp > .5*max(max(tmp));
%imagesc(2.^TF_vec00,flipud(2.^SF_vec00),flipud(tmp3'))
MaxSF00 = s.x(4);
MaxTF00 = s.x(5);

indSF50 = find(k2b00>.5*s.Maxfit00 & x00(:,1)>MaxSF00);
indTF50 = find(k2b00>.5*s.Maxfit00 & x00(:,2)>MaxTF00);
SFhicut_50 = NaN;
TFhicut_50 = NaN;
if ~isempty(indSF50)
    SFhicut_50 = max(x00(indSF50,1));
end
if ~isempty(indTF50)
    TFhicut_50 = max(x00(indTF50,2));
end
s.SFhicut_50 = SFhicut_50;
s.TFhicut_50 = TFhicut_50;

%repeat for 10% cutoff

indSF10 = find(k2b00>.1*s.Maxfit00 & x00(:,1)>MaxSF00);
indTF10 = find(k2b00>.1*s.Maxfit00 & x00(:,2)>MaxTF00);
SFhicut_10 = NaN;
TFhicut_10 = NaN;
if ~isempty(indSF10)
    SFhicut_10 = max(x00(indSF10,1));
end
if ~isempty(indTF10)
    TFhicut_10 = max(x00(indTF10,2));
end
s.SFhicut_10 = SFhicut_10;
s.TFhicut_10 = TFhicut_10;




if PLOTIT_FIT == 1
    figure;
    subplot(2,2,1);
    imagesc(s.data); colormap('gray'); axis image; colorbar;
    subplot(2,2,2);
    imagesc(s.k2b); colormap('gray'); axis image;colorbar;
    subplot(2,2,4);
    imagesc(s.res); colormap('gray'); axis image;colorbar;
end

if SAVEALLDATA == 0
    s = rmfield(s,'res');
    s = rmfield(s,'k2b');
    s = rmfield(s,'k2');
    s = rmfield(s,'data');
end



%this.gab1 = s;



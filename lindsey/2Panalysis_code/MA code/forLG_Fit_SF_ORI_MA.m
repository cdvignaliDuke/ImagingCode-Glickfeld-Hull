%Fitting ori x SF: 

%PLAN:
%0) plot raw data, ori vs SF
%1) average 2 dirns, then find best ori (avg'ed over SF).
%2) extract SF params using DOG, using best ori
%3) separately, find best SF (at pref ori), then fit ori tuning curve at that
%best SF
%4) repeat steps 2-3 100 times after sampling with replacement, get bootstrap
%estimates for SF, SF_width, ori and ori_width

s = struct;
s.data = data;


%first, collapse across directions: 
data2 = (data(1:4,:) + data(5:8,:))./2;

%m = @(pars,dummy) Gauss2D_ellipseMA(pars,grid2.sfsf,grid2.tftf);


%step 1: 
%

%find approx best ori: 
data3 = mean(data2,2);
ind_bestori = find(data3 == max(data3));
ind_bestori = ind_bestori(1);
s.ind_bestori = ind_bestori;
s.data3 = data3;

%AK CODE: 
%SF_Y = response
%SF_X = spatial frequencies in cpd
%SF_vec0,ORI_vec0
SF_X = SF_vec0';
if SF_X(1) == 0
    SF_X(1) = SF_X(2)/2; %if SF_X(2) == .02, then fullfield is set to .01
end

SF_Y = data2(ind_bestori,:); %NOTE: could add data at fullfield later.. 
%CHANGING THIS: to improve SNR, using all oris: 
%SF_Y = mean(data2(:,:),1); %NOTE: could add data at fullfield later.. 


%NOTE: for now, make all negative values = 0 (could do this better,
%potentially, by subtracting the blank stimulus.. 

s.SF_X = SF_X;
s.SF_Y = SF_Y;



options.Display = 'off';
[b,r,J,COVB,mse] = nlinfit(SF_X,SF_Y,@DOGfn,[0.15 0.15 6 6],options);
[ypred,delta] = nlpredci(@DOGfn,SF_X,b,r,'covar',COVB);
             
m = @(pars,X) DOGfn(pars,X);
s.m = func2str(m);
s.fit = b; %temp(ind);
%s.x = s.fit.x2;
%s.k2 = m(s.fit.x2,x);
s.k2b = ypred; %reshape(s.k2,m2,n2);
s.res = s.SF_Y-s.k2b';


%%%%%%%%%
%now find SF peak and bandwidth from the DOG fit
%SF_X2 = 2.^[-7:.1:-1]';

%only fit between min and max:
%SF_X2 = 2.^[-7:.1:-1]';
Min = log2(SF_X(1));
Max = log2(SF_X(end));

SF_X2 = 2.^[Min:.1:Max]; 

Y = DOGfn(b,SF_X2);
[i,j] = max(Y);
SFmax = SF_X2(j(1));

if SFmax >= min(SF_X) & SFmax <= max(SF_X)
    s.SFmax = SFmax;
elseif SFmax < min(SF_X)
    s.SFmax = min(SF_X) - .0001;
elseif SFmax > max(SF_X);
    s.SFmax = max(SF_X)+.0001;
end

%find bandwith
ind1 = 1:j(1);
ind_lower = max(find(Y(ind1)<(.5*max(Y))));

ind1 = j(1):length(Y);
ind_upper = min(find(Y(ind1)<(.5*max(Y)))) + j(1) - 1;

if ~isempty(ind_lower) & ~isempty(ind_upper)
    SFband = log2(Y(ind_upper)) - log2(Y(ind_lower));
    s.SFband = SFband;
else
    s.SFband = NaN;
end
    

%now choose the peak SF and fit ori tuning curve (and then estimate
%direction tuning..
%ORI_X = ORI_vec0';
[i,j] = max(data2(ind_bestori,:));
%ORI_Y = data2(:,j(1)); %NOTE: could add data at fullfield later.. 

%now  compute classical ori estimates

i = 1;
OriStatKO(i).ori_ratio_change = data2(:,j(1));
OriStatKO(i).dir_ratio_change = data(:,j(1));
% classical analysis
[OriStatKO(i).best_dir,...
    OriStatKO(i).null_dir,...
    OriStatKO(i).R_best_dir,...
    OriStatKO(i).R_null_dir,...
    OriStatKO(i).R_min_dir,...
    OriStatKO(i).DI]...
    = dir_indexKO(OriStatKO(i).dir_ratio_change);
[OriStatKO(i).best_ori,...
    OriStatKO(i).null_ori,...
    OriStatKO(i).R_best_ori,...
    OriStatKO(i).R_null_ori,...
    OriStatKO(i).R_min_ori,...
    OriStatKO(i).OI]...
    = dir_indexKO(OriStatKO(i).ori_ratio_change);

% vector average analysis

[OriStatKO(i).dir_vector_angle,...
    OriStatKO(i).dir_vector_mag,...
    OriStatKO(i).dir_vector_tune]...
    = vector_average(OriStatKO(i).dir_ratio_change);

[OriStatKO(i).ori_vector_angle,...
    OriStatKO(i).ori_vector_mag,...
    OriStatKO(i).ori_vector_tune]...
    = vector_average(OriStatKO(i).ori_ratio_change);

OriStatKO(i).ori_vector_angle = OriStatKO(i).ori_vector_angle /2;

s.OriStatKO = OriStatKO;


%s.x: SFband SFmax oritune oriangle dirtune dirangle

%old way: 
%s.x = [s.SFband s.SFmax ...
%    OriStatKO(i).ori_vector_tune OriStatKO(i).ori_vector_angle ...
%    OriStatKO(i).dir_vector_tune OriStatKO(i).dir_vector_angle];

OI2 = (OriStatKO(i).R_best_ori - OriStatKO(i).R_null_ori)./(OriStatKO(i).R_best_ori + OriStatKO(i).R_null_ori);
DI2 = (OriStatKO(i).R_best_dir - OriStatKO(i).R_null_dir)./(OriStatKO(i).R_best_dir + OriStatKO(i).R_null_dir);

%s.x = [s.SFband s.SFmax ...
%    OriStatKO(i).OI OriStatKO(i).ori_vector_angle ...
%    OriStatKO(i).DI OriStatKO(i).dir_vector_angle];
s.x = [s.SFband s.SFmax ...
    OI2 OriStatKO(i).ori_vector_angle ...
    DI2 OriStatKO(i).dir_vector_angle];

return



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


%A sigma_SF sigma_TF sf0 tf0 xi
s.lb =  [.001 .25 .25  (min(x(:,1))-1) (min(x(:,2))-1) -2];
s.ub =  [1000 3   3    (max(x(:,1))+1) (max(x(:,2))+1) 2];
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
s.k2 = m(s.fit.x2,x);
s.k2b = reshape(s.k2,m2,n2);
s.res = s.data-s.k2b;

if PLOTIT_FIT == 1
    figure(gcf); clf
    subplot(2,2,1);
    imagesc(s.data'); colormap('gray'); axis image;
    subplot(2,2,2);
    imagesc(s.k2b'); colormap('gray'); axis image;
    subplot(2,2,4);
    imagesc(s.res'); colormap('gray'); axis image;
end

%this.gab1 = s;


return;

x = -30:30;
[xx,yy]=meshgrid(x,x);

sf = .1;
ori = 45;
pars =         [.001 .5/6 .005 -pi/2 -30 -30 -pi/2];
pars =         [1000  4/6 .1    pi/2  30  30  pi/2];
pars = [15   2/6 .02   4*pi/4   0   0     0];
%pars = [1 1/sf sf (0) 0 0 (ori/180*pi)];
h = gabor(pars,xx,yy);
figure;imagesc(h);


%pars_start = [15   2/6 .02   pi/4   0   0     0];
pars_start = [15   1/6 .02   0   0   0     0];

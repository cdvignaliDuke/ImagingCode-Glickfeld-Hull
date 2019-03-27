%VB Gabor fitting, modified to fit 2D elliptical

%params VB sends in: sta.kmean

s = struct;
s.data = data;
s.orig = reshape(b',size(b,1)*size(b,2),1);

m = @(pars,azel) Gauss2D_ellipseMA(pars,azel);
s.m = func2str(m);

y = s.data(:);

CM = zeros(1,2);
CM_data = zeros(1,2);
CM(1) = sum(x_plot(:,1).*s.orig)/sum(s.orig);
CM(2) = sum(x_plot(:,2).*s.orig)/sum(s.orig);
CM_data = CM;


Mask_NaN = isnan(y);
ind_NaN = find(isnan(y));
ind_noNaN = find(Mask_NaN == 0);
if ~isempty(ind_NaN)
    x_plot(ind_NaN,:) = [];
    y(ind_NaN) = [];
end



%A sigma_Az sigma_El Az0 El0 xi
s.lb =  [.001 1 1  (min(x_plot(:,1))) (min(x_plot(:,2))) 0];
s.ub =  [3  100  100    (max(x_plot(:,1))) (max(x_plot(:,2))) 0];



Nsamps = 2; %changed from 2 to 3
dbin = [s.ub(1) - s.lb(1); ...
    s.ub(2) - s.lb(2); ...
    s.ub(3) - s.lb(3); ...
    s.ub(4) - s.lb(4); ...
    s.ub(5) - s.lb(5)]./ (Nsamps-1);
%s.ub(6) - s.lb(6)] ./ (Nsamps-1);

sigma_El_vec = s.lb(3)+dbin(3)/2:dbin(3):s.ub(3);
sigma_Az_vec = s.lb(2)+dbin(2)/2:dbin(2):s.ub(2);
El_vec = s.lb(5)+dbin(5)/2:dbin(5):s.ub(5);
Az_vec = s.lb(4)+dbin(4)/2:dbin(4):s.ub(4);
%xi_vec = s.lb(6)+dbin(6)/2:dbin(6):s.ub(6);


s.x0 =  [max(max(s.data)) sigma_Az_vec sigma_El_vec Az_vec El_vec 0]; %xi_vec(ixi)];
options = optimset('Display', 'off');
[x2,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]  = lsqcurvefit(m,s.x0,x_plot,y,s.lb,s.ub,options);

s.SStot = sum(sum((s.data-mean(s.data)).^2));
s.Rsq = 1 - Resnorm/s.SStot;

s.x = x2;
s.k2 = zeros(m2*n2,1);
s.k2(ind_noNaN) = m(x2,x_plot);
s.k2_plot = m(x2,x_plot);

s.k2_plot_oversamp0 = m(s.x,xyhigh);
s.k2_plot_oversamp = reshape(s.k2_plot_oversamp0,xNperfreq,yNperfreq);

s.k2b = reshape(s.k2,m2,n2);
s.k2b_plot = reshape(s.k2_plot,m2,n2);

s.res = s.data-s.k2b;
s.Maxfit = max(max(s.k2b));
s.Maxdata = max(max(s.data));

s.fval = FVAL;

%find highcut for SF and TF
%x(:,1) = log2(grid2.sfsf(:));
%x(:,2) = log2(grid2.tftf(:));
%    TF_vec0 = ind_TFuse(:,2);
%    SF_vec0 = ind_SFuse(:,2);
x00 = zeros(size(grid2.ElEl00,1)*size(grid2.ElEl00,2),2);
x00(:,1) = (grid2.AzAz00(:));
x00(:,2) = (grid2.ElEl00(:));
k2b00 = m(s.x,x00);
%now find 10% and 50%
s.Maxfit00 = max(max(k2b00));

%tmp2 = max(tmp,[],1); %max over TFs
%tmp3 = tmp > .5*max(max(tmp));
%imagesc(2.^TF_vec00,flipud(2.^SF_vec00),flipud(tmp3'))
MaxAz00 = s.x(4);
MaxEl00 = s.x(5);


indAz50 = find(k2b00>.5*s.Maxfit00 & x00(:,1)>MaxAz00);
indEl50 = find(k2b00>.5*s.Maxfit00 & x00(:,2)>MaxEl00);
Elhicut_50 = NaN;
Azhicut_50 = NaN;
if ~isempty(indEl50)
    Elhicut_50 = max(x00(indEl50,2));
end
if ~isempty(indAz50)
    Azhicut_50 = max(x00(indAz50,1));
end
s.Elhicut_50 = Elhicut_50;
s.Azhicut_50 = Azhicut_50;

%repeat for 10% cutoff
indAz10 = find(k2b00>.1*s.Maxfit00 & x00(:,1)>MaxAz00);
indEl10 = find(k2b00>.1*s.Maxfit00 & x00(:,2)>MaxEl00);

Elhicut_10 = NaN;
Azhicut_10 = NaN;
if ~isempty(indEl10)
    Elhicut_10 = max(x00(indEl10,2));
end
if ~isempty(indAz10)
    Azhicut_10 = max(x00(indAz10,1));
end
s.Elhicut_10 = Elhicut_10;
s.Azhicut_10 = Azhicut_10;



if PLOTIT_FIT == 1
    if start ==65
        set(gcf, 'Position', [0 0 800 1000]);
        UN = username;
        if UN(1:3) == 'kev'
            fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RFfits' num2str(ifig) '.pdf']);
        else
            fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_RFfits' num2str(ifig) '.pdf']);
        end
        print(fn_out,'-dpdf')
        figure;
        ifig = 1+ifig;
        start = 1;
    end
    h = subplot(8,8,start);
    
    MIN = 0;
    MAX = max(max(data,[],1),[],2);
    dF = num2str(MAX*100);
    imagesc(s.data); colormap('gray'); axis image;
    if h_all(1,iCell)
        sig_str = ' **';
    else
        sig_str = [];
    end
    title([num2str(iCell) sig_str]);
    set(h,'XTick',[1:size(uvar.Az,1)])
    set(h,'YTick',[1:size(uvar.El,1)])
    set(h,'YTickLabel',flipud(uvar.El));set(h,'XTickLabel',(uvar.Az))
    axis off
    h = subplot(8,8, 1+start);
    imagesc(s.k2b_plot); colormap('gray'); axis image;
    set(h,'XTick',[1:size(uvar.Az,1)])
    set(h,'YTick',[1:size(uvar.El,1)])
    set(h,'YTickLabel',flipud(uvar.El));set(h,'XTickLabel',(uvar.Az))
    axis off
    h = subplot(8,8, 2+start);
    imagesc(s.k2_plot_oversamp); colormap('gray'); axis image;
    axis off
    h = subplot(8,8,3+start);
    imagesc(s.res); colormap('gray'); axis image;
    set(h,'XTick',[1:size(uvar.Az,1)])
    set(h,'YTick',[1:size(uvar.El,1)])
    set(h,'YTickLabel',flipud(uvar.El));set(h,'XTickLabel',(uvar.Az))
    axis off
    clim([0 MAX])
    start = start+4;
end

if SAVEALLDATA == 0
    s = rmfield(s,'res');
    s = rmfield(s,'k2b');
    s = rmfield(s,'k2');
    s = rmfield(s,'data');
end



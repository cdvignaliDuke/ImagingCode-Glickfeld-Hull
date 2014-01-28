m = @(pars,angle) Gauss_dir_LG(pars,angle);
s.x0 =  [min(data(:,ind_bestSFTF)) OriStatKO.R_best_dir OriStatKO.R_null_dir s.y(OriStatKO.best_dir).*pi./180 pi./2];
s.ub = [3 3 3 2*pi 2*pi];
s.lb = [0 0 0 0 pi/4];
options = optimset('Display', 'off');
[x2,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]  = lsqcurvefit(m,s.x0,s.y ,data(:,ind_bestSFTF)',s.lb,s.ub,options);
temp(index) = cell2struct({x2,Resnorm},names,2);

function [best_ori_fit, tuning_width, A1, k1, A, resnorm, residual, exitflag, output]...
    = estimate_ori_tuningKO (ydata, DoInterp, method, printflag)

% estimate tuning paramteters by fitting with single von Mises function
% ydata: average signal value for each direction
% DoInterp: if 0, no interpolation. if 1, ydata is interpolated with double sampling.
% method: interpolation method. e.g. spilne, cubic, linear...
%
% tuning_width: orientation tuning width from curve fitting. 
% A1: amplitude of the peak for orientation (ratio change).
% k1: inverse (not exact) of orientation tuning width.
% best_ori_fit: preferred orientation obtained form curve fitting. 
% resnorm: residual of curve fitting for orientation.
%
%   Kenichi Ohki 09/16/04
%



if nargin < 2
    DoInterp = 0;
end
if nargin < 3
    method = 'linear';
end
if nargin < 4
    printflag = 0;
end

nstim_per_run = length(ydata);

xdata = [0:nstim_per_run-1] * 180 / nstim_per_run;

[ymax, best_ori] = max(ydata);
ymin = min(ydata);
ymin=0;
% first guess
A0(1) = ymax - ymin;        %A1
A0(2) = 1;                %k1
A0(3) = (best_ori -1)*180/nstim_per_run; %phi1

options = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on');

% interpolate data to double directions

ydatai = [ydata(nstim_per_run),ydata,ydata(1:2)];
ydatai = interp1([-1:nstim_per_run+1],ydatai,[-1:0.5:nstim_per_run+1],method);
ydatai = ydatai(3:nstim_per_run*2+2);
xdatai= [0:nstim_per_run*2-1] * 90 / nstim_per_run;

if DoInterp == 0;
    [A,resnorm,residual,exitflag,output] = lsqcurvefit(@vonMises1, A0, xdata, ydata-ymin, [], [], options);
else
    [A,resnorm,residual,exitflag,output] = lsqcurvefit(@vonMises1, A0, xdatai, ydatai-ymin, [], [], options);
end

if printflag == 1;
    figure;
    x=[0:180];
    plot(x, VonMises1(A, x)+ymin);
    hold on;
    plot([0:nstim_per_run]*180/nstim_per_run, [ydata(1:nstim_per_run),ydata(1)], 'r');
    hold off;
    print ('-dpsc', '-append', 'fitting.ps');

end

A1=A(1);
k1=A(2);
best_ori_fit=mod(A(3),180);

  
if k1>=log(2)/2
    tuning_width = acos(log(0.5)/k1+1)/pi*180/2;    
else
    tuning_width= 90;
end

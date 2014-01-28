function [best_dir_fit, null_dir_fit, R_best_dir_fit, R_null_dir_fit, DI_fit, DS,...
        tuning_width, A1, A2, k1, k2, phi2, A, resnorm,...
        residual, exitflag, output] = estimate_dir_tuningKO (ydata, DoInterp, method, printflag)

% estimate tuning paramteters by fitting with double von Mises function
% ydata: average signal value for each direction
% DoInterp: if 0, no interpolation. if 1, ydata is interpolated with double sampling.
% method: interpolation method. e.g. spilne, cubic, linear...
%
% best_dir_fit: preferred direction obtained from curve fitting.
% DIfit: direction index from curve fitting
%               (1-R(best_dir_fit+180)/R(best_dir_fit)).
% DS: Swindale's direction selectivity index.
% tuning_width: direction tuning width from curve fitting.
% A1: amplitude of the larger peak for direction (ratio change).
% A2: amplitude of the smaller peak for direction (ratio change).
% k1: inverse (not exact) of direction tuning width for the larger peak. 
% k2: inverse (not exact) of direction tuning width for the smaller peak. 
% phi2: direction of the second peak.
% resnorm: residual of curve fitting for direction.
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

xdata = [0:nstim_per_run-1] * 360 / nstim_per_run;

[ymax, best_dir] = max(ydata);
ymin = min(ydata);
ymin=0;
null_dir = mod((best_dir + nstim_per_run/2 - 1), nstim_per_run)+1;

% first guess
A0(1) = ymax - ymin;        %A1
A0(2) = ydata(null_dir)-ymin;    %A2
A0(3) = 2.0;                %k1
A0(4) = 2.0;                %k2
A0(5) = (best_dir -1)*360/nstim_per_run; %phi1
A0(6) = (null_dir -1)*360/nstim_per_run; %phi2

options = optimset('LargeScale', 'off', 'LevenbergMarquardt', 'on');

% interpolate data to double directions

ydatai = [ydata(nstim_per_run),ydata,ydata(1:2)];
ydatai = interp1([-1:nstim_per_run+1],ydatai,[-1:0.5:nstim_per_run+1],method);
ydatai = ydatai(3:nstim_per_run*2+2);
xdatai= [0:nstim_per_run*2-1] * 180 / nstim_per_run;

if DoInterp == 0;
    [A,resnorm,residual,exitflag,output] = lsqcurvefit(@vonMises2, A0, xdata, ydata-ymin, [], [], options);
else
    [A,resnorm,residual,exitflag,output] = lsqcurvefit(@vonMises2, A0, xdatai, ydatai-ymin, [], [], options);
end

if printflag == 1;
    figure;
    x=[0:360];
    plot(x, VonMises2(A, x)+ymin);
    hold on;
    plot([0:nstim_per_run]*360/nstim_per_run, [ydata(1:nstim_per_run),ydata(1)], 'r');
    hold off;
    print ('-dpsc', '-append', 'fitting.ps');
 
end

A(5)= mod(A(5),360);
A(6)= mod(A(6),360);

   
if A(1) >= A(2);
    A1=A(1);
    A2=A(2);
    k1=A(3);
    k2=A(4);
    best_dir_fit=A(5);
    phi2=A(6);
else
    A1=A(2);
    A2=A(1);
    k1=A(4);
    k2=A(3);
    best_dir_fit=A(6);
    phi2=A(5);
end

null_dir_fit = mod(best_dir_fit + 180, 360);
R_best_dir_fit = VonMises2(A, best_dir_fit)+ymin;
R_null_dir_fit = VonMises2(A, null_dir_fit)+ymin;

DI_fit = 1-R_null_dir_fit/R_best_dir_fit;
DS = (A1-A2)/(A1+A2);

if k1>=log(2)/2
    tuning_width = acos(log(0.5)/k1+1)/pi*180;    
else
    tuning_width = 180;
end

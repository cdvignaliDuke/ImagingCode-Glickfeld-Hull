function [best_dir_interp, null_dir_interp, R_best_dir_interp, R_null_dir_interp, DI_interp] = interp_dir_tuningKO (ydata, method, printflag)

% interpolate direction tuning curve
% ydata: average signal value for each direction
% method: interpolation method. e.g. spilne, cubic, linear...
% best_dir_interp: preferred direction from interpolation.
% DI_interp: direction index from interpolation
%               (1-R(best_dir_interp)/R(best_dir_interp)).
%
%   Kenichi Ohki 09/22/04
%


if nargin < 2
    method = 'linear';
end
if nargin < 3
    printflag = 0;
end

nstim_per_run = length(ydata);


% interpolate data to 360 points

if method == 'fourier'
    ydatai = fourier_interp (ydata);
else
    ydatai = [ydata(nstim_per_run),ydata,ydata(1:2)];
    ydatai = interp1([-1:nstim_per_run+1],ydatai,[-1:nstim_per_run/360:nstim_per_run+1],method);
    ydatai = ydatai(360/nstim_per_run+1:360/nstim_per_run+360);
end

[R_best_dir_interp, best_dir_interp] = max(ydatai);

best_dir_interp = best_dir_interp -1;
null_dir_interp = mod(best_dir_interp+180, 360);

R_null_dir_interp = ydatai(null_dir_interp+1);

DI_interp = 1-R_null_dir_interp/R_best_dir_interp;

if printflag == 1;
    figure;
    x=[0:359];
    plot(x, ydatai);
    hold on;
    plot([0:nstim_per_run]*360/nstim_per_run, [ydata(1:nstim_per_run),ydata(1)], 'r');
    hold off;
    print ('-dpsc', '-append', 'interp.ps');
 
end

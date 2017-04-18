function [yfit, b] = slr(x, y_obs, win)
% get simple linear regression on x and y_obs within time-window, win. x
% and y are both column vectors
% yfit is fit line, b is slope and intersept

X = x(win);
Y = y_obs(win);

b = polyfit(X,Y,1);
yfit = polyval(b,X,1);

end

function plot_cone_input (resp,i,j,name)

% Sept. 2008  Kenichi Ohki 

x=resp(i,:);
y=resp(j,:);
plot(x,y,'.')
hold on
title(name)
A=min([x,y]);
B=max([x,y]);
plot([A,B],[0,0]);
plot([0,0],[A,B]);
plot([A,B],[A,B]);
axis ([A,B,A,B])
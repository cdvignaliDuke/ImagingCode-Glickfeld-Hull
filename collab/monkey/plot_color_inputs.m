
function plot_color_inputs(resp)

% resp: 10 colors x N cells matrix
% it can be dF/F or dF/F/contrast or dF
% Sept. 2008  Kenichi Ohki 

LonMon=1;
Mon=2;
MonLoff=3;
Loff=4;
LoffMoff=5;
Moff=6;
LonMoff=7;
Lon=8;
Son=9;
Soff=10;

figure

subplot(2,3,1);
plot_cone_input(resp,Son,Soff,'S');

subplot(2,3,2);
plot_cone_input(resp,Mon,Moff,'M');

subplot(2,3,3);
plot_cone_input(resp,Lon,Loff,'L');

subplot(2,3,4);
plot_cone_input(resp,LonMoff,MonLoff,'L-M');

subplot(2,3,5);
plot_cone_input(resp,LonMon,LoffMoff,'L+M');





function plot_cone_contrast (resp, filter1, filter2, filter3)

% resp: 10 colors x N cells matrix
% it can be dF/F or dF/F/contrast or dF
% Sept. 2008  Kenichi Ohki 

Ncells=size(resp,2)


if nargin<2
    filter1=[1:Ncells];
end

if nargin<3
    filter2=[];
end

if nargin<4
    filter3=[];
end



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


for i=1:Ncells
    S(i)=resp(Son,i)-resp(Soff,i);
    M(i)=resp(Mon,i)-resp(Moff,i);
    L(i)=resp(Lon,i)-resp(Loff,i);
    Sum(i)=abs(S(i))+abs(M(i))+abs(L(i));
    Sn(i)=S(i)./Sum(i);
    Mn(i)=M(i)./Sum(i);
    Ln(i)=L(i)./Sum(i);
end

Ncells = 452;

figure
plot(Mn(filter1),Ln(filter1),'.r');
hold on;
plot(Mn(filter2),Ln(filter2),'.g');
plot(Mn(filter3),Ln(filter3),'.b');

plot([-1,0],[0,1],'k');
plot([-1,0],[0,-1],'k');
plot([1,0],[0,1],'k');
plot([1,0],[0,-1],'k');
plot([-1,1],[0,0],'k');
plot([0,0],[-1,1],'k');

axis([-1,1,-1,1])

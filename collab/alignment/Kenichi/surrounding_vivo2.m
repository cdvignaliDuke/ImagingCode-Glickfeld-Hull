% Find points in centroidsv that  surround a chosen point B. 
% input: the chosen points B=[x y z]
% output: points (by findcell3D) surrounding B; TCv are the coordinates.

% label surrounding cells within a radius r
prompt = {'x:','y:','z:','r'};
dlg_title = 'label surrounding cells';
answer=inputdlg(prompt,dlg_title);
A=[str2num(answer{1}),str2num(answer{2}),str2num(answer{3})];
r=str2num(answer{4});

Sv =zeros(5,100);
j=0;
for i=[1:size(centroidsv,2)];
    if abs(centroidsv(1,i)-A(1))<r && abs(centroidsv(2,i)-A(2))<r && abs(centroidsv(3,i)-A(3))<r
        j=j+1;  
        Sv(:,j)=cat(1, centroidsv(1:3,i),1,i); 
    end;
end;

% Overlape the result on fixedcentroidsv.fig.
hold on;
scatter3(Sv(1,1:j),Sv(2,1:j),Sv(3,1:j),50,'b','filled');

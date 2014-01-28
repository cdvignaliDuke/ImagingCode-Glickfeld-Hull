% Find points in tformed-centroidsf that  surround a chosen point A. 
% input: the chosen points A=[x y z]
% output: points (by findcell3D) surrounding A; TCf are the coordinates.


% label surrounding cells within a radius r
prompt = {'x:','y:','z:','r'};
dlg_title = 'label surrounding cells';
answer=inputdlg(prompt,dlg_title);
A=[str2num(answer{1}),str2num(answer{2}),str2num(answer{3})];
r=str2num(answer{4});


Sftc =zeros(5,100);
j=0;
for i=[1:size(tformedCPDed,2)];
    if abs(tformedCPDed(1,i)-A(1))<r && abs(tformedCPDed(2,i)-A(2))<r && abs(tformedCPDed(3,i)-A(3))<r
        j=j+1;
        Sftc(:,j)=cat(1, tformedCPDed(1:3,i),1,i); 
    end
end

% Overlape the result on fixedtformed.fig.
hold on;
scatter3(Sftc(1,1:j),Sftc(2,1:j),Sftc(3,1:j),50,'m','filled');

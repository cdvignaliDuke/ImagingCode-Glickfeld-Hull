function pos=lblMatrix2pos(lblMatrix, resolution)
%convers nD label matix to position of each label
% position is matrix index
% if resolution is provided then position is in um in X, Y, Z order
% and top left pixel has (0,0) coordinate

lMax=max(lblMatrix(:));

dims=length(size(lblMatrix));



tmpstr=[];

for n=1:dims
    if n>1
        tmpstr=[tmpstr ','];
    end
    tmpstr=[tmpstr 'A' num2str(n) ' '];
end
pos=zeros(lMax,dims);
for lN=1:lMax
    ind=find(lblMatrix==lN);
    eval(['[' tmpstr ']=ind2sub(size(lblMatrix),ind);']);
    for d=1:dims
        tmpstr2=['A' num2str(d)'];
        eval(['pos(lN,d)=mean(' tmpstr2 ');']);
    end
end

if nargin>1 
    pos_tmp=pos-1;
    pos(:,1)=pos_tmp(:,2)*resolution(1);
    pos(:,2)=pos_tmp(:,1)*resolution(2);
    if dims>2
        pos(:,3)=pos_tmp(:,3)*resolution(3);     
    end
end
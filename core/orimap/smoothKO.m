function out = smoothKO(in)

dim=size(in);

out=zeros(dim(1),1);

out(1)=in(1)/2+in(2)/4+in(dim(1))/4;

for i=2:dim(1)-1
    out(i)=in(i)/2+in(i-1)/4+in(i+1)/4;
end

out(dim(1))=in(dim(1))/2+in(dim(1)-1)/4+in(1)/4;



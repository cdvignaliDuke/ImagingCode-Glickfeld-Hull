function [out1 out2 out3]=transform(coeffs,in1,in2,in3)

if nargin~=2 && nargin~=4
    error('transform - number of input parameters should be 2 or 4');
end

if nargin==2
    sz=size(in1);
    out1=zeros(sz);
    out2=[];
    out3=[];
    inx=in1(:,1);
    iny=in1(:,2);
    inz=in1(:,3);
else
    inx=in1;
    iny=in2;
    inz=in3;
end

outx=coeffs(1)+coeffs(2)*inx+coeffs(3)*iny+coeffs(4)*inz;
outy=coeffs(5)+coeffs(6)*inx+coeffs(7)*iny+coeffs(8)*inz;
outz=coeffs(9)+coeffs(10)*inx+coeffs(11)*iny+coeffs(12)*inz;


if nargin==2
    out1(:,1)=outx;
    out1(:,2)=outy;
    out1(:,3)=outz;
else
    out1=outx;
    out2=outy;
    out3=outz;
end
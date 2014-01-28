function hueKO=hueHSV2KO(hueHSV)

% convert hue values (0-1) in HSV color space to hue vaules in KO color
% space
%
% Kenichi Ohki  2008.11.4

temp=hueHSV;
temp(find(temp<0))=0;
temp(find(temp>1))=1;
a=find(temp<=1/3);
b=find(temp>1/3);
temp(a)=temp(a)*3/2;
temp(b)=temp(b)*3/4+1/4;

hueKO=temp;




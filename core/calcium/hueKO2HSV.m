function hueHSV=hueKO2HSV(hueKO)

% convert hue values (0-1) in KO color space to hue vaules in HSV color
% space
%
% Kenichi Ohki  2008.11.4

% temp=hueKO;
% temp(find(temp<0))=0;
% temp(find(temp>1))=1;
% a=find(temp<=0.5);
% b=find(temp>0.5);
% temp(a)=temp(a)*2/3;
% temp(b)=temp(b)*4/3-1/3;
% 
% hueHSV=temp;

hueHSV=hueKO;
hueHSV(hueHSV<0)=0;
hueHSV(hueHSV>1)=1;
hueHSV=hueHSV*2/3;
hueHSV(hueHSV>(1/3))=hueHSV(hueHSV>(1/3))*2-1/3;


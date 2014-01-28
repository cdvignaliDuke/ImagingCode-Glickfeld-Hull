function combo_coeffs=calculateComboTransf(coeffs1,coeffs2)
combo_coeffs=[];

points=rand(20,3);

[outx1 outy1 outz1]=transform(coeffs1,points(:,1),points(:,2),points(:,3));
[outx2 outy2 outz2]=transform(coeffs2,outx1,outy1,outz1);

Xr = [ones(size(points,1),1) points(:,1) points(:,2) points(:,3)];
warning off all
combo_coeffs(1:4)=regress(outx2,Xr);
combo_coeffs(5:8)=regress(outy2,Xr);
combo_coeffs(9:12)=regress(outz2,Xr);
warning on all


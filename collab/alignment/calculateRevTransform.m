function rev_coeffs=calculateRevTransform(coeffs)

points=rand(20,3);

% calculate inverse transformation
[outx outy outz]=transform(coeffs,points(:,1),points(:,2),points(:,3));
%calculate transformation using regression
Xr = [ones(size(points,1),1) outx outy outz];
warning off all
rev_coeffs(1:4)=regress(points(:,1),Xr);
rev_coeffs(5:8)=regress(points(:,2),Xr);
rev_coeffs(9:12)=regress(points(:,3),Xr);
warning on all
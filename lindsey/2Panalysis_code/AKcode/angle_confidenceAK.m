function CI = angle_confidenceAK (a,alpha)
    
[ori_boot_th,r] = cart2pol(mean(sind(a)),mean(cosd(a)));
ori_boot_angle = 90-atand(ori_boot_th);
n=length(a);
R=n*r;
X=chi2inv((1-alpha),1);

if ((r<=0.9)&(r>(X/(2*n))^(1/2)))
    CI=acosd((1/R)*((2*n*(2*R^2-n*X))/(4*n-X))^(1/2));
else
    if r>0.9
        CI=acosd((1/R)*(n^2-(n^2-R^2)*exp(X/n))^(1/2));
    else
        CI=NaN;
    end
end
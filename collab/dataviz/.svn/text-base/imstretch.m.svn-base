function img3 = imstretch(img,tol,gam)

if nargin < 2
    tol = [0.001 0.999];
end;

if nargin < 3
    gam = 1;
end

img2 = imscale(img);
img3 = imadjust(img2,stretchlim(img2,tol),[],gam);
imshow(img3)

return;


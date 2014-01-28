function playstim(imgs,n,fps)
%PLAYSTIM
% PLAYSTIM(IMGS,N,FPS)

if nargin < 2
    n = 1;
end

if nargin < 3
    fps = 30;
end

[ny,nx,nframes]=size(imgs);

clear mov;
for ind = 1:nframes
    pic = double(imgs(:,:,ind))+1;
    mov(ind)= im2frame(pic,gray(256));
end

imshow(imgs(:,:,1));
tic;
movie(mov,n,fps);
t = toc;

fprintf(1,'Actual refresh rate = %2.1f fps\n',n*nframes/t);

return;
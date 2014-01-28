function [dir_ratio, ori_ratio] = calc_ratio_changeb (dir_dF_sm, ori_dF_sm, bases_sm)

% calculate dF/F images to each direction and orientation.

dim=size(dir_dF_sm);
ndir=dim(3);
nori=ndir/2;

dir_ratio=zeros(dim(1), dim(2), ndir);
ori_ratio=zeros(dim(1), dim(2), nori);


for i=1:ndir
    dir_ratio(:,:,i)=dir_dF_sm(:,:,i)./bases_sm(:,:,i);
end

for i=1:nori
    ori_ratio(:,:,i)=(dir_ratio(:,:,i)+dir_ratio(:,:,i+nori))/2;
end


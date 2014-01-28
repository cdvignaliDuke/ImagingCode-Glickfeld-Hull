function [dir_dF, dir_dF_sm, ori_dF, ori_dF_sm] = calc_dF (dir_F, base, sp_filter)

% calculate dF images (and their smoothed images) to each direction and orientation.
% if sp_filter is omitted, dir_dF_sm = ori_dF_sm = 0.

dim=size(dir_F);
ndir=dim(3); %default=3
nori=ndir/2;
dir_dF=zeros(dim(1), dim(2), ndir);
dir_dF_sm=zeros(dim(1), dim(2), ndir);
ori_dF=zeros(dim(1), dim(2), nori);
ori_dF_sm=zeros(dim(1), dim(2), nori);


for i=1:ndir
    dir_dF(:,:,i)=dir_F(:,:,i)-base;
    if (nargin >= 3)
        dir_dF_sm(:,:,i)=filter2(sp_filter, dir_dF(:,:,i));
    else
        dir_dF_sm(:,:,i)=0;
    end
end

for i=1:nori
    ori_dF(:,:,i)=(dir_dF(:,:,i)+dir_dF(:,:,i+nori))/2;
    if (nargin >= 3)
        ori_dF_sm(:,:,i)=(dir_dF_sm(:,:,i)+dir_dF_sm(:,:,i+nori))/2;
    else
        ori_dF_sm(:,:,i)=0;
    end
end


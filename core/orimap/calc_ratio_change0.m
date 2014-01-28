function dir_ratio = calc_ratio_change0 (dir_dF_sm, base_sm)

% calculate dF/F images to each direction.

dim=size(dir_dF_sm);
ndir=dim(3);

dir_ratio=zeros(dim(1), dim(2), ndir);


for i=1:ndir
    dir_ratio(:,:,i)=dir_dF_sm(:,:,i)./base_sm;
end

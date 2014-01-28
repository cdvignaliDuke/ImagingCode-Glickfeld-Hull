function ave = stack_ave (stack)

[xDim, yDim, tDim]=size(stack);

ave=zeros(xDim,yDim);    

for t=1:tDim
    ave=ave+stack(:,:,t);
end

ave = ave ./ tDim;


function stack_out = contrastnorm_3DMA(stack, sigma, IS_rectanglestack, NoDataMask);
    %contrast normalize:
    %stack = 3D volume
    %sigma indicates degree of convolution, in pixels, 
    %IS_rectanglestack: 0 -> square image, 1 -> 1x2 rectange image, so set
    %    the rest of thesquare image to the mean of the valid part. 
    %NoDataMask: optional input to avoid smoothing pixels that are zero due
    %    to corregistration shifting.. 
    % -MA Feb09
    
    stack3_vol_avgB = double(stack);
    
    if nargin < 4
        ISMASK = 0;
    else
        ISMASK = 1;
    end
    if nargin < 3
        IS_rectanglestack = 0; %default square
    end
    
    %sigma = 10; %pixels: 
    
    X1 = [-2*sigma:(2*sigma)];
    
    [X,Y] = ndgrid(X1);
    Gauss2D = exp(-(X.^2 + Y.^2)./sigma^2);
    Gauss2D = Gauss2D./sum(sum(sum(Gauss2D)));

    Lx = size(stack3_vol_avgB,1);
    Ly = size(stack3_vol_avgB,2);
    Lz = size(stack3_vol_avgB,3);
    Nconv = length(X1)/2-1/2;
    tmp0 = zeros(Lx+2*Nconv, Ly + 2*Nconv);
    
   stack3_vol_avgC= zeros(size(stack3_vol_avgB));
    for count = 1:size(stack3_vol_avgB,3);
        tmpMask = zeros(Lx,Ly);
        if IS_rectanglestack == 0
            tmpMask(:) = 1;
        else
            tmpMask(1:Lx/2,1:Ly) = 1;
        end
        if ISMASK == 1
            ind = find(squeeze(NoDataMask(:,:,count)));
            tmpMask(ind) = 0;
        end
        ind_USE = find(tmpMask);
        ind_OMIT = find(tmpMask==0);
        tmp = squeeze(stack3_vol_avgB(:,:,count));


        tmp0(:) = mean(mean(tmp(ind_USE)));
        tmp_new = tmp;
        tmp_new(ind_OMIT) = mean(mean(tmp(ind_USE)));
        if IS_rectanglestack == 0
            tmp0((Nconv+1):(Nconv+Lx),(Nconv+1):(Nconv+Ly)) = tmp_new;
         else
            tmp0((Nconv+1):(Nconv+Lx/2),(Nconv+1):(Nconv+Ly)) = tmp_new(1:Lx/2,1:(Ly));
        end
        
        tmp2 = tmp ./ conv2(tmp0,Gauss2D,'valid');;
        stack3_vol_avgC(:,:,count) = tmp2;
    end
 stack_out = stack3_vol_avgC;


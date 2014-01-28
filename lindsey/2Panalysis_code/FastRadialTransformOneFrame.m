function [pos ra ma mi debugImg]=FastRadialTransformOneFrame(imgin, radii1,radii2, posGuess,motion_range,contrast_tsh,alpha, kappa,subsampl,smooth,squeeze_coeff)
% by SY model after Loy & Zrlinsky
% 2011/04/04
%imgin - input image
%radii1 - range of reflection radii to test
%radii2 - range of pupl radii to test
%posGuess - initial estimation of reflection and pupl position [reflX reflY piplX puplY]
%motion_range - possible deviation of reflection and pupl from initial position estimation
%contrast_tsh - treshold for image gradient (small gradient is ignored)
%alpha - method parameter
%kappa - method parameter
%subsampl - image subsample 
%smooth - smoothing parameter for postprocessing
%squeeze_coeff - compensation for image elongation, for horiztontaly elongated images squeeze_coeff<1


if nargin<7 || isempty(alpha)
   alpha=10;
end

if nargin<8 || isempty(kappa)
    kappa=10;
end

if nargin<9 || isempty(subsampl)
    subsampl=1;
end

if nargin<10 || isempty(smooth)
    smooth=1;
end

if nargin<11 || isempty(squeeze_coeff)
    squeeze_coeff=1;
end

if subsampl>1
    small_img=imgin(1:subsampl:end,1:subsampl:end);
    radii1=radii1/subsampl;
    radii2=radii2/subsampl;
    posGuess=posGuess/subsampl;
    motion_range=motion_range/subsampl;
    smooth=smooth/subsampl;
else
    small_img=imgin;
end



sz=size(small_img);
[gradX,gradY] = gradient(small_img); % use built-in Matlab function instead of Sobel
mag=sqrt(gradX.^2 +gradY.^2);
gradX=gradX./mag;
gradY=gradY./mag;

[X,Y] = meshgrid(1:sz(2),1:sz(1));


S1=zeros(sz(1),sz(2));
S2=zeros(sz(1),sz(2));


me=median(small_img(:));

% cut number of pixels to consider (two groups)
ind2=find(mag(:)>contrast_tsh & small_img(:)<me);
ind1=find(mag(:)>contrast_tsh & small_img(:)>me);



maRa1=max(radii1);
maRa2=max(radii2);
Omax=zeros(length(radii1),1);


% filters
sp_filter=fspecial('gaussian', round(smooth*2)+1, smooth);
f_sz=size(sp_filter);
[meshX, meshY]=meshgrid(-(f_sz(1)-1)/2:(f_sz(1)-1)/2, -(f_sz(2)-1)/2:(f_sz(2)-1)/2);
sp_filter((meshX.^2+meshY.^2)>(f_sz(1)/2).^2)=0;
sp_filter=sp_filter/sum(sp_filter(:)); %renormalize filter (sum of all filter pixels

sp_filter_2=fspecial('gaussian', round(smooth*2*2)+1, smooth*2);
f_sz=size(sp_filter_2);
[meshX, meshY]=meshgrid(-(f_sz(1)-1)/2:(f_sz(1)-1)/2, -(f_sz(2)-1)/2:(f_sz(2)-1)/2);
sp_filter_2((meshX.^2+meshY.^2)>(f_sz(1)/2).^2)=0;
sp_filter_2=sp_filter/sum(sp_filter_2(:)); %renormalize filter (sum of all filter pixels

% search for reflection
for n=1:length(radii1)
    r=radii1(n);
    O=zeros(sz(1),sz(2));
    
    rX=round(r*gradX);
    rY=round(r*gradY*squeeze_coeff);
    posX=X+rX;
    posY=Y+rY;
    
    if isempty(posGuess)
        for t=1:length(ind1)
            if posY(ind1(t))>maRa1+2 && posY(ind1(t))<sz(1)-maRa1-2 && posX(ind1(t))>maRa1+2 && posX(ind1(t))<sz(2)-maRa1-2 ...
                    && posY(ind1(t))>0 && posY(ind1(t))<sz(1)+1 && posX(ind1(t))>0 && posX(ind1(t))<sz(2)+1
                O(posY(ind1(t)),posX(ind1(t)))=O(posY(ind1(t)),posX(ind1(t)))+1;
            end
        end
    else
        for t=1:length(ind1)
            if posY(ind1(t))>posGuess(2)-motion_range(1) && posY(ind1(t))<posGuess(2)+motion_range(1) ...
                    && posX(ind1(t))>posGuess(1)-motion_range(1) && posX(ind1(t))<posGuess(1)+motion_range(1) ...
                    && posY(ind1(t))>0 && posY(ind1(t))<sz(1)+1 && posX(ind1(t))>0 && posX(ind1(t))<sz(2)+1
                O(posY(ind1(t)),posX(ind1(t)))=O(posY(ind1(t)),posX(ind1(t)))+1;
            end
        end
    end
    Omax(n)=max(O(:));
    O=min(O,kappa);
    F1=(O./kappa).^alpha;
    S1=S1+F1;
end

Omin=zeros(length(radii2),1);
for n=1:length(radii2)
    r=radii2(n);
    O=zeros(sz(1),sz(2));
    
    rX=round(r*gradX);
    rY=round(r*gradY*squeeze_coeff);
    negX=X-rX;
    negY=Y-rY;
    
    if isempty(posGuess)
        for t=1:length(ind2)
            if negY(ind2(t))>maRa2+2 && negY(ind2(t))<sz(1)-maRa2-2 && negX(ind2(t))>maRa2+2 && negX(ind2(t))<sz(2)-maRa2-2 ...
                    && posY(ind2(t))>0 && posY(ind2(t))<sz(1)+1 && posX(ind2(t))>0 && posX(ind2(t))<sz(2)+1
                O(negY(ind2(t)),negX(ind2(t)))=O(negY(ind2(t)),negX(ind2(t)))-1;
            end
        end
    else
        for t=1:length(ind2)
            if negY(ind2(t))>posGuess(4)-motion_range(2) && negY(ind2(t))<posGuess(4)+motion_range(2) ...
                    && negX(ind2(t))>posGuess(3)-motion_range(2) && negX(ind2(t))<posGuess(3)+motion_range(2) ...
                    && negY(ind2(t))>0 && negY(ind2(t))<sz(1)+1 && negX(ind2(t))>0 && negX(ind2(t))<sz(2)+1
                O(negY(ind2(t)),negX(ind2(t)))=O(negY(ind2(t)),negX(ind2(t)))-1;
            end
        end
    end
%    O=conv2(double(O),sp_filter,'same'); %improves accuracy of radius calculation
    Omin(n)=min(O(:));
    O=max(O,round(-1.2*kappa));
    %F2=sign(O).*(abs(O)./kappa).^alpha;
    F2=-(-O./kappa).^alpha;
    S2=S2+F2;
end



% puple - using center of mass
[smi smiind]=min(S2(:));
mind=(S2(:)<smi*0.3);
miX=sum(X(:).*mind(:).*S2(:))/sum(S2(mind));
miY=sum(Y(:).*mind(:).*S2(:))/sum(S2(mind));

% reflection - using smoothed max position
 S=conv2(double(S1),sp_filter,'same');
 [ma maind]=max(S(:));
 [maY maX]=ind2sub(sz,maind);

if smi==0
    miX=posGuess(3);
    miY=posGuess(4);
    mi=0;
else
    % calculate "goodness of fit for center  position only
    S2ext=zeros(sz+smooth*2);
    S2ext(smooth+1:end-smooth,smooth+1:end-smooth)=S2;
    S4=S2ext(round(miY)-smooth+smooth:round(miY)+smooth+smooth,round(miX)-smooth+smooth:round(miX)+smooth+smooth).*sp_filter_2;        
    mi=sum(S4(:));
end

pos=([maX maY miX miY]-1)*subsampl+1;




[rma rmaind]=max(Omax(:));
[rmi0 rmiind]=min(Omin(:));

ra(1)=radii1(rmaind);

rmind=Omin(:)<(rmi0*0.7);
ra(2)=sum(radii2(:).*rmind.*Omin(:))/sum(Omin(rmind));

if isnan(ra(2))
    ra(2)=radii2(round(mean(1:length(radii2))));
end

ra=ra*subsampl;


debugImg=[];
% im=mand*2+mind;
% im=reshape(im,sz);
%debugImg=uint8((S1(1:2/subsampl:end,1:2/subsampl:end))/sma*255); %assuming subsample is 1 or 2


function [r_m_out_tr g_m_out_tr b_m_out_tr]=readTIFandtransform(fname,transform_coeff,limits,resolution_in,resolution_out)
%[r_m_out_tr g_m_out_tr
%b_m_out_tr]=readTIFandtransform('E:\users\Sergey\Matlab\Lindsey\091215\invivostacks\091215_invivo_RGB_um_rotate-45_fliphor_reslicetop_86to185.tif',transform_coeff,[512 512 49],[0.34 0.34 1.02],[0.27 0.27 2]);
%trnasforms tif stack according to transormantion coeffs

r_m_out_tr=[];
g_m_out_tr=[];
b_m_out_tr=[];

inv_transform_coeff=calculateRevTransform(transform_coeff);


% read tif file info
info = imfinfo(fname,'tif');

if strcmp(info(1).ColorType, 'grayscale')
    colorImg=0;
else
    colorImg=1;
end


W=info(1).Width;
H=info(1).Height;
nFrames=length(info);
tmp_im= imread(fname,1);
data_class=class(tmp_im(1));

if nargin<3 || isempty(limits)
%    limits=[512 512 256];

Xcorners=([1 1 H H 1 1 H H]-1)*resolution_in(1);
Ycorners=([1 W 1 W 1 W 1 W]-1)*resolution_in(2);
Zcorners=([1 1 1 1 nFrames nFrames nFrames nFrames]-1)*resolution_in(3);
[Xout,Yout,Zout]=transform(transform_coeff,Xcorners,Ycorners,Zcorners);


if round(min(Xout)/resolution_out(1))<0 || round(min(Yout)/resolution_out(2))<0 || round(min(Zout)/resolution_out(3))<0
    display ('Stack is shifted');
    limits=[round(max(Xout)/resolution_out(1))+1 round(max(Yout)/resolution_out(2))+1 round(max(Zout)/resolution_out(3))+1];
    limitsmin=[round(min(Xout)/resolution_out(1))+1 round(min(Yout)/resolution_out(2))+1 round(min(Zout)/resolution_out(3))+1];
else
    limits=[round(max(Xout)/resolution_out(1))+1 round(max(Yout)/resolution_out(2))+1 round(max(Zout)/resolution_out(3))+1];
    limitsmin=[1 1 1];
end
else
    
    limitsmin=[1 1 1];
    
end

if nargin<4
    resolution_in=[1 1 1];
end

if nargin<5
    resolution_out=[1 1 1];
end


%preallocate color stacks
% Red_im=zeros(H,W,nFrames,data_class);
% if colorImg
%     Green_im=zeros(H,W,nFrames,data_class);
%     Blue_im=zeros(H,W,nFrames,data_class);
% end


Red_im=zeros(H,W,nFrames);
if colorImg
    Green_im=zeros(H,W,nFrames);
    Blue_im=zeros(H,W,nFrames);
end

%read frames, separete color chennels
for frame=1:nFrames
    tmp_im= imread(fname,frame);
    Red_im(:,:,frame)=double(tmp_im(:,:,1));
    if colorImg
        Green_im(:,:,frame)=double(tmp_im(:,:,2));
        Blue_im(:,:,frame)=double(tmp_im(:,:,3));
    end
end

tic

maxXout=limits(1);
maxYout=limits(2);
maxZout=limits(3);

if isempty(limitsmin)
    minXout=1;
    minYout=1;
    minZout=1;
else
    minXout=limitsmin(1);
    minYout=limitsmin(2);
    minZout=limitsmin(3);
end

%[mshXout0,mshYout0,mshZout0] = meshgrid(0:resolution_out(1):(maxXout-1)*resolution_out(1),0:resolution_out(2):(maxYout-1)*resolution_out(2),0:resolution_out(3):(maxZout-1)*resolution_out(3));


[mshXout0,mshYout0,mshZout0] = meshgrid(((minXout:maxXout)-1)*resolution_out(1),((minYout:maxYout)-1)*resolution_out(2),((minZout:maxZout)-1)*resolution_out(3));
[mshXout,mshYout,mshZout]=transform(inv_transform_coeff,mshXout0,mshYout0,mshZout0);
[mshXin,mshYin,mshZin] = meshgrid(0:resolution_in(1):(W-1)*resolution_in(1),0:resolution_in(2):(H-1)*resolution_in(2),0:resolution_in(3):(nFrames-1)*resolution_in(3));



clear mshXout0;
clear mshYout0;
clear mshZout0;


r_m_out_tr = interp3(mshXin,mshYin,mshZin,Red_im,mshXout,mshYout,mshZout,'linear');
clear Red_im;

if colorImg
    g_m_out_tr = interp3(mshXin,mshYin,mshZin,Green_im,mshXout,mshYout,mshZout,'linear');
    clear Green_im;
    b_m_out_tr = interp3(mshXin,mshYin,mshZin,Blue_im,mshXout,mshYout,mshZout,'linear');
    clear Blue_im;
end
toc
% 
% %preallocate color stacks
% Red_im=zeros(H,W,nFrames,data_class);
% if colorImg
%     Green_im=zeros(H,W,nFrames,data_class);
%     Blue_im=zeros(H,W,nFrames,data_class);
% end
% 
% 
% %read frames, separete color chennels
% for frame=1:nFrames
%     tmp_im= imread(fname,frame);
%     Red_im(:,:,frame)=tmp_im(:,:,1);
%     if colorImg
%         Green_im(:,:,frame)=tmp_im(:,:,2);
%         Blue_im(:,:,frame)=tmp_im(:,:,3);
%     end
% end



% r_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Red_im),mshXout,mshYout,mshZout,'linear');
% if colorImg
%     g_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Green_im),mshXout,mshYout,mshZout,'linear');
%     b_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Blue_im),mshXout,mshYout,mshZout,'linear');
% end




out_size=size(r_m_out_tr);
fname_out=strrep(fname,'.','_transform.');

clear Red_im;
clear Green_im;
clear Blue_im;

for frame=1:out_size(3)
    if colorImg
        warning off all
        RGB_frame=zeros(out_size(1),out_size(2),3,data_class);
        RGB_frame(:,:,1)=r_m_out_tr(:,:,frame);
        RGB_frame(:,:,2)=g_m_out_tr(:,:,frame);
        RGB_frame(:,:,3)=b_m_out_tr(:,:,frame);
        warning on all
    else
        eval(['RGB_frame=' data_class '(r_m_out_tr(:,:,frame));']);
    end

    if frame==1
        imwrite(RGB_frame,fname_out,'WriteMode','overwrite','Compression','none');
    else
        imwrite(RGB_frame,fname_out,'WriteMode','append','Compression','none');
    end

end


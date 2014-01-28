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

if nargin<3
    limits=[512 512 256];
end

if nargin<4
    resolution_in=[1 1 1];
end

if nargin<5
    resolution_out=[1 1 1];
end



maxXout=limits(1);
maxYout=limits(2);
maxZout=limits(3);



[mshXout0,mshYout0,mshZout0] = meshgrid(0:resolution_out(1):(maxXout-1)*resolution_out(1),0:resolution_out(2):(maxYout-1)*resolution_out(2),0:resolution_out(3):(maxZout-1)*resolution_out(3));

[mshXin,mshYin,mshZin] = meshgrid(0:resolution_in(1):(W-1)*resolution_in(1),0:resolution_in(2):(H-1)*resolution_in(2),0:resolution_in(3):(nFrames-1)*resolution_in(3));

[mshXout,mshYout,mshZout]=transform(inv_transform_coeff,mshXout0,mshYout0,mshZout0);


clear mshXout0;
clear mshYout0;
clear mshZout0;

%preallocate color stacks
Red_im=zeros(H,W,nFrames,data_class);
if colorImg
    Green_im=zeros(H,W,nFrames,data_class);
    Blue_im=zeros(H,W,nFrames,data_class);
end


%read frames, separete color chennels
for frame=1:nFrames
    tmp_im= imread(fname,frame);
    Red_im(:,:,frame)=tmp_im(:,:,1);
    if colorImg
        Green_im(:,:,frame)=tmp_im(:,:,2);
        Blue_im(:,:,frame)=tmp_im(:,:,3);
    end
end

r_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Red_im),mshXout,mshYout,mshZout,'linear');
if colorImg
    g_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Green_im),mshXout,mshYout,mshZout,'linear');
    b_m_out_tr = interp3(mshXin,mshYin,mshZin,double(Blue_im),mshXout,mshYout,mshZout,'linear');
end


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


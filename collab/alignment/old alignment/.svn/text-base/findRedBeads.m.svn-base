function coord_B=findRedBeads(fname,resolution_in)
%coord_B=finddRedBids('E:\users\Sergey\Matlab\Lindsey\091215\invivostacks\091215_invivo_RGB_um_rotate-45_fliphor_reslicetop_86to185.tif',transform_coeff,[0.34 0.34 1.02]);


coord_B=[];
if nargin<2
    resolution_in=[1 1 1];
end

tsh=0.8;
max_v=round(120/(resolution_in(1)*resolution_in(2)*resolution_in(3)));

% read tif file info
info = imfinfo(fname,'tif');

W=info(1).Width;
H=info(1).Height;
nFrames=length(info);
tmp_im= imread(fname,1);
data_class=class(tmp_im(1));


%preallocate stacks
Red_im=zeros(H,W,nFrames);

%read only red chennels
for frame=1:nFrames
    tmp_im= imread(fname,frame);
    Red_im(:,:,frame)=tmp_im(:,:,1);
end


%figure;
m1=max(Red_im,[],1);
m2=max(m1,[],2);
%plot(squeeze(m2));

mm=zeros(length(m2),5);
mm(:,1)=m2;
mm(1:end-1,2)=m2(2:end);
mm(1:end-2,3)=m2(3:end);
mm(3:end,4)=m2(1:end-2);
mm(2:end,5)=m2(1:end-1);

mmm=max(mm,[],2);
ma=max(Red_im(:));
mmm=max(mmm,ma*0.3);
ma3 = shiftdim(mmm,-2);
ma4=repmat(ma3,[H W 1]);


best_tsh=0;
bestLN=0;
for tsh=0.3:0.02:1
    im_bw=Red_im>tsh*ma4;
    BW_L = bwareaopen(im_bw,8,6);
    L = bwlabeln(BW_L,6);
    LN=max(L(:));
    if LN>=bestLN
        bestLN=LN;
        best_tsh=tsh;
    end
end

best_tsh
im_bw=Red_im>best_tsh*ma4;
BW_L = bwareaopen(im_bw,8,6);


BW_L_2 = bwareaopen(BW_L,max_v,6);
BW_S= BW_L & ~BW_L_2;

Red_im_2=Red_im.*BW_L_2;

best_tsh_2=0;
bestLN=0;
for tsh=best_tsh:0.01:1
    im_bw_2=Red_im_2>tsh*ma4;
    BW_L_3 = bwareaopen(im_bw_2,8,6);
    L = bwlabeln(BW_L_3,6);
    LN=max(L(:));
    if LN>=bestLN
        bestLN=LN;
        best_tsh_2=tsh;
    end
end

best_tsh_2
im_bw_2=Red_im_2>best_tsh_2*ma4;
BW_L_3 = bwareaopen(im_bw_2,8,6);


BW_L=BW_S | BW_L_3;


L = bwlabeln(BW_L,6);


fname_out=strrep(fname,'.','_bw.');
for frame=1:nFrames
    if frame==1
        imwrite(BW_L(:,:,frame),fname_out,'WriteMode','overwrite','Compression','none');
    else
        imwrite(BW_L(:,:,frame),fname_out,'WriteMode','append','Compression','none');
    end
end

% BW_L = bwareaopen(im_bw,max_v,6);
% BW_S= im_bw & ~BW_L;
% 
% 
% fname_out=strrep(fname,'.','_bw_l.');
% for frame=1:nFrames
%     if frame==1
%         imwrite(BW_L(:,:,frame),fname_out,'WriteMode','overwrite','Compression','none');
%     else
%         imwrite(BW_L(:,:,frame),fname_out,'WriteMode','append','Compression','none');
%     end
% end
% 
% 
% fname_out=strrep(fname,'.','_bw_s.');
% for frame=1:nFrames
%     if frame==1
%         imwrite(BW_S(:,:,frame),fname_out,'WriteMode','overwrite','Compression','none');
%     else
%         imwrite(BW_S(:,:,frame),fname_out,'WriteMode','append','Compression','none');
%     end
% end


%L = bwlabeln(BW_S,26);
%L = bwlabeln(im_bw,6);

% fname_out=strrep(fname,'.','_bw_le.');
% for frame=1:nFrames
%     if frame==1
%         imwrite(L(:,:,frame),fname_out,'WriteMode','overwrite','Compression','none');
%     else
%         imwrite(L(:,:,frame),fname_out,'WriteMode','append','Compression','none');
%     end
% end


stat = regionprops(L, 'Centroid');
coord_B=cat(1, stat.Centroid)-1;

coord_B(:,1)=coord_B(:,1)*resolution_in(1);
coord_B(:,2)=coord_B(:,2)*resolution_in(2);
coord_B(:,3)=coord_B(:,3)*resolution_in(3);

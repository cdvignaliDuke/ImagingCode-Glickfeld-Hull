function out=Shift4D(stack,target,Nplanes,is_rect,reverse_ch,Nupsamp,Nfromedge,Nvolbin)

IS_rectanglestack = is_rect; %1 for e.g. 128x256 images embedded in a 256x256 matrix, 0 for 256x256 or other square image.. 
REVERSE_CH1CH2_VERT = reverse_ch; % if recording blue and green channels, then green is on bottom, need to move it to top for consistency
SKIPLASTSCANLINE = 1; %skips processing of last scan line on fastscan axis due to large increases in intensity
USETARG = 1;

Lx = size(stack,1);
Ly = size(stack,2);
Lz = size(stack,3);
Nvol = floor(Lz./Nplanes);


if SKIPLASTSCANLINE == 1;
    if IS_rectanglestack == 1; %1 for e.g. 128x256 images embedded in a 256x256 matrix, 0 for 256x256 or other square image.. 
        stack(Lx/2,:,:) = 0;
    elseif IS_rectanglestack == 0; 
        stack(Lx,:,:) = 0;
    end
end

%now reshape into 4D stack: 
%stack4D = zeros(Lx,Ly,Lz, Nvol);
stack4D_USE = double(reshape(stack(:,:,1:Nplanes*Nvol),Lx,Ly,Nplanes,Nvol));
clear('stack');
if REVERSE_CH1CH2_VERT == 1 & IS_rectanglestack
    stack4D_USE([(Lx/2 + [1:Lx/2]) 1:Lx/2],:,:,:) = stack4D_USE;
end


%for simplicity, remove the 1st plane due to flyback, for now
stack4D_USE(:,:,1,:) = [];
target(:,:,1)=[];
Lz = size(stack4D_USE,3);



%stack4D_USE =  squeeze(stack4D(Xuse,Yuse,:,1:100));

%upsample the Zaxis? 
Nupsamp_ZONLY = 1;
if Nupsamp_ZONLY > 1
    Z00a = 1:Lz;
    
    [X00 Y00 Z00] = meshgrid(1:Lx,1:Ly,Z00a);
    
    X0a = 1:Lx;
    Y0a = 1:Ly;
    Z0a = 1:1/Nupsamp_ZONLY:Lz;

    Lx2 = length(X0a);
    Ly2 = length(Y0a);
    Lz2 = length(Z0a);
    [X0 Y0 Z0] = meshgrid(1:Lx, 1:Ly, Z0a);
    stack4D_USE2 = zeros(Lx2,Ly2,Lz2, Nvol);
    tmp = squeeze(stack4D_USE(:,:,:,1));
    tmp2 = reshape(tmp,Lx*Ly,Lz)';
    tmp3 = interp1(Z00a,tmp2,Z0a);
    
    Vol_targ2 = zeros(Lx2,Ly2,Lz2);
    
    FASTVERSION = 1; % '1' -> use interp1, '3' -> use interp3
    if FASTVERSION == 1
        display(['upsampling Z dimension by ', num2str(Nupsamp_ZONLY)]); tic
        for count = 1:Nvol
%            fprintf([num2str(count),'\n']);
            fprintf([num2str(count),'_']);
            tmp(:) = squeeze(stack4D_USE(:,:,:,count));
            tmp2(:) = reshape(tmp,Lx*Ly,Lz)';
            tmp3(:) = interp1(Z00a,tmp2,Z0a);
            stack4D_USE2(:,:,:,count) = reshape(tmp3',Lx2,Ly2,Lz2); %interp3(X00, Y00, Z00,stack4D_USE(:,:,:,count),X0, Y0, Z0);
        end
        tmp(:) = single(target);
        tmp2(:) = reshape(tmp,Lx*Ly,Lz)';
        tmp3(:) = interp1(Z00a,tmp2,Z0a);
        Vol_targ2 = reshape(tmp3',Lx2,Ly2,Lz2);
        display(['done upsampling Z dimension by ', num2str(Nupsamp_ZONLY)]); toc
    else
        display(['upsampling Z dimension by ', num2str(Nupsamp_ZONLY)]); tic
        for count = 1:Nvol
            fprintf([num2str(count),'\n']);
            stack4D_USE2(:,:,:,count) = interp3(X00, Y00, Z00,stack4D_USE(:,:,:,count),X0, Y0, Z0);
        end
        Vol_targ2(:) = interp3(X00, Y00, Z00,target,X0, Y0, Z0);
        display(['done upsampling Z dimension by ', num2str(Nupsamp_ZONLY)]); toc
    end
    target = single(Vol_targ2);
    clear('stack4D_USE');
    stack4D_USE = stack4D_USE2;
    clear('stack4D_USE2');
end

target=single(target);

%now take avg: 
%stack_avg0(:,:,1) = [];
%stack_avg0 = squeeze(mean(stack4D_USE,4));

%now for corregistration: 

%use a smaller part of XYplane for correg: 
%Xuse = floor(Lx/4):floor(Lx*3/4);
%Yuse = floor(Ly/4):floor(Ly*3/4);

%stack4D_USE =  squeeze(stack4D(Xuse,Yuse,:,1:100));
%stack4D_USE =  squeeze(stack4D_USE(:,:,:,:));
%clear('stack4D');

%%stack_avg00 = squeeze(mean(stack4D_USE,4));
stack_avg00 = target;
%Lx2 = length(Xuse);
%Ly2 = length(Yuse);
ZOOMIN = 0;

%choose subset of image to use for estimating shifts: 
if ZOOMIN == 1
    if IS_rectanglestack == 1
%        WIN_IM = floor([[1 (Lx/2 - 1)]; [Ly/4 Ly*3/4]; [1 size(stack_avg00,3)]]); %part of volume to use, in X, Y, Z dim
        WIN_IM = floor([[Lx/8 (3*Lx/8)]; [Ly/4 Ly*3/4]; [1 size(stack_avg00,3)]]); %part of volume to use, in X, Y, Z dim
    else
        WIN_IM = floor([[Lx/4 Lx*3/4]; [Ly/4 Ly*3/4]; [1 size(stack_avg00,3)]]); %part of volume to use, in X, Y, Z dim
    end
elseif ZOOMIN == 0
    if IS_rectanglestack == 1
        WIN_IM = [[1 (Lx/2 - 1)]; [1 Ly]; [1 size(stack_avg00,3)]]; %part of volume to use, in X, Y, Z dim
    else
        WIN_IM = [[1 (Lx - 1)]; [1 Ly]; [1 size(stack_avg00,3)]]; %part of volume to use, in X, Y, Z dim
    end
end
%NOTE: in above, we're excluding the last line in the X axis from analysis
%because of changes in scan intensity in the last line.. 

%Nupsamp = [1]; %no subpixel = 1, subpixel -> integer>1
%Nfromedge = [4 4 4]; %distance in pixels for potential correg movement, in X, Y, Z ([4 4 2] means +-3, +-3, +-1 pixels)
%Nvolbin = 16; %SET TO 16 for AK expts; for correg step, avg all data from Nvolbin volumes, then correg and apply that shift to each of the Nvolbin volumes 

%imagesc(stack4D_USE(:,:,6,2))

ind_vol_USE = 1:Nvol;
%ind_vol_USE = 1:size(stack4D_USE,4);
%[stack4D_reg, ShiftMat] = reg3D_xcorrMA(stack4D_USE(:,:,:,ind_vol_USE), stack_avg00, WIN_IM, Nupsamp, Nfromedge, Nvolbin);

Nvolbin=1;
[ShiftMat ShiftMat_subp] = reg3D_xcorrMA_shiftsonly_subp(stack4D_USE(:,:,:,ind_vol_USE), target, WIN_IM, Nupsamp, Nfromedge, Nvolbin);

%realign volumes, upsampling in Z
%Nupsamp_XYZ = [1 1 2]; %number of times to upsample registered data in x, y, z
%[stack4D_reg, Mask_postreg] = reg3D_xcorrMA_realignonly_subp(stack4D_USE(:,:,:,ind_vol_USE), ShiftMat_subp,Nupsamp_XYZ);
%stack4D_reg(find(stack4D_reg == -1)) = 0;

%unique(stack4D_reg)
%out=cell(1,2);
%out{1,1}= uint16(stack4D_reg);
%out{1,2} = ShiftMat_subp;
out=ShiftMat_subp;
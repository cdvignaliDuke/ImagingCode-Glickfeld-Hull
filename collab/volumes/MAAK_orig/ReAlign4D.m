function out=ReAlign4D(stack,ShiftMat,Nplanes,is_rect,reverse_ch,Nupsamp_XYZ)

IS_rectanglestack = is_rect; %1 for e.g. 128x256 images embedded in a 256x256 matrix, 0 for 256x256 or other square image.. 
REVERSE_CH1CH2_VERT = reverse_ch; % if recording blue and green channels, then green is on bottom, need to move it to top for consistency
SKIPLASTSCANLINE = 1; %skips processing of last scan line on fastscan axis due to large increases in intensity

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

ind_vol_USE = 1:Nvol;

%realign volumes, upsampling in Z
%Nupsamp_XYZ = [1 1 2]; %number of times to upsample registered data in x, y, z
[stack4D_reg, Mask_postreg] = reg3D_xcorrMA_realignonly_subp(stack4D_USE(:,:,:,ind_vol_USE),ShiftMat,Nupsamp_XYZ);
stack4D_reg(find(stack4D_reg == -1)) = 0;

out= uint16(stack4D_reg);
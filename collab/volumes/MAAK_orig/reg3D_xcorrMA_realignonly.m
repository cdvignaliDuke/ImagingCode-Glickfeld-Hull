
%reg3D_xcorrMA.m
%%now do 3D corregistration
function [vol_avg_shift_down] = reg3D_xcorrMA_realignonly( vol_avg_USE2, ShiftMat, Nupsamp, Nfromedge, Nvolbin );
%Nfromedge should be all even numbers


%WIN_IM = [1 128; 1 size(vol_avg,2); 1 (size(vol_avg,3)-1)];
%Nupsamp = 4;
%Nfromedge = [4 4 3]; 



if nargin < 3
    Nupsamp = 1;
end
if nargin < 4
    Nfromedge = [4 4 4];
end
if nargin < 5
    Nvolbin = 1;
end

Nvol0 = size(vol_avg_USE2,4);
if Nvolbin == 1
%    vol_avg = vol_avg0;
%    clear('vol_avg0');
    Nvol = Nvol0;
else
    Nvol = floor(Nvol0./Nvolbin);
end


Lx2 = size(vol_avg_USE2,1);
Ly2 = size(vol_avg_USE2,2);
Lz2 = size(vol_avg_USE2,3);


NfromedgeX = Nfromedge(1);
NfromedgeY = Nfromedge(2);
NfromedgeZ = Nfromedge(3); %starting pixel in Z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now realign using estimated shifts

  display('realigning volumes'); tic

%Lx2 = size(vol_avg,1);
%Ly2 = size(vol_avg,2);
%Lz2 = size(vol_avg,3);

%NfromedgeX2 = 1*Nupsamp;
%NfromedgeY2 = 1*Nupsamp;
%NfromedgeZ2 = 1*Nupsamp;

%find out max shift

Lx3 = Lx2 + 2*(NfromedgeX-1);
Ly3 = Ly2 + 2*(NfromedgeY-1);
Lz3 = Lz2 + 2*(NfromedgeZ-1);

Im_enlarged = zeros(Lx3,Ly3,Lz3);

ind_X_shift0 = (NfromedgeX):(Lx2+NfromedgeX-1);
ind_Y_shift0 = (NfromedgeY):(Ly2+NfromedgeY-1);
ind_Z_shift0 = (NfromedgeZ):(Lz2+NfromedgeZ-1);

%Lx2_shift = length(ind_X_shift0);
%Ly2_shift = length(ind_Y_shift0);
%Lz2_shift = length(ind_Z_shift0);

%vol_avg_shift = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol0);

vol_avg_shift = zeros(Lx2,Ly2,Lz2,Nvol0);



%vol_avg_shift_XYonly = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol);
%vol_avg_noshift = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol);
%vol_avg_shift_AVG = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol);
%vol_avg_shift_XYonly_AVG = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol);
%vol_avg_noshift_AVG = zeros(Lx2_shift,Ly2_shift,Lz2_shift,Nvol);
%vol_avg_USE2_AVG = squeeze(mean(vol_avg_USE2,4));

for count = 1:Nvol
    if Nvolbin == 1
        ind_X_shift = ind_X_shift0 + ShiftMat(count,1);
        ind_Y_shift = ind_Y_shift0 + ShiftMat(count,2);
        ind_Z_shift = ind_Z_shift0 + ShiftMat(count,3);

        vol_avg_shift0 = zeros(Lx3, Ly3, Lz3);
        vol_avg_shift0(ind_X_shift, ind_Y_shift, ind_Z_shift) = ...
            vol_avg_USE2(:,:,:,count);
        vol_avg_shift(:,:,:,count) = vol_avg_shift0(ind_X_shift0, ...
            ind_Y_shift0, ...
            ind_Z_shift0);
        
    elseif Nvolbin > 1
        ind_X_shift = ind_X_shift0 + ShiftMat(count*Nvolbin,1);
        ind_Y_shift = ind_Y_shift0 + ShiftMat(count*Nvolbin,2);
        ind_Z_shift = ind_Z_shift0 + ShiftMat(count*Nvolbin,3);

            vol_avg_shift0 = zeros(Lx3, Ly3, Lz3, Nvolbin);
        vol_avg_shift0(ind_X_shift, ind_Y_shift, ind_Z_shift,:) = ...
            vol_avg_USE2(:,:,:,(count-1)*Nvolbin + [1:Nvolbin]);
        vol_avg_shift(:,:,:,(count-1)*Nvolbin + [1:Nvolbin]) = vol_avg_shift0(ind_X_shift0, ...
            ind_Y_shift0, ind_Z_shift0,:);
        
    end

       
 % vol_avg_shift_XYonly(:,:,:,count) = vol_avg_USE2(ind_X_shift, ...
%					    ind_Y_shift, ...
%					    ind_Z_shift0,count);
 % vol_avg_noshift(:,:,:,count) = vol_avg_USE2(ind_X_shift0, ...
%					    ind_Y_shift0, ...
%					    ind_Z_shift0,count);
   %NOW do the same using the avg volume
%  vol_avg_shift_AVG(:,:,:,count) = vol_avg_USE2_AVG(ind_X_shift, ...
%					    ind_Y_shift, ...
%					    ind_Z_shift);
%  vol_avg_shift_XYonly_AVG(:,:,:,count) = vol_avg_USE2_AVG(ind_X_shift, ...
%					    ind_Y_shift, ...
%					    ind_Z_shift0);
%  vol_avg_noshift_AVG(:,:,:,count) = vol_avg_USE2_AVG(ind_X_shift0, ...
%					    ind_Y_shift0, ...
%					    ind_Z_shift0);
   
                    
end



  %downsample back down to original X and Y size: 

  Ndownsamp = Nupsamp;
  
  ind_X_down = 1:Ndownsamp:Lx2;
  ind_Y_down = 1:Ndownsamp:Ly2;
  ind_Z_down = 1:Ndownsamp:Lz2;

 % vol_avg_shift_down = vol_avg_shift(ind_X_down,ind_Y_down,ind_Z_down,:);
  if Nupsamp > 1
      vol_avg_shift_down = zeros(size(vol_avg));
      vol_avg_shift_down(:,:,:,:) = vol_avg_shift(ind_X_down,ind_Y_down,ind_Z_down,:);
  else
      vol_avg_shift_down = vol_avg_shift;
  end
  
  display('done realigning volumes'); toc

      clear('vol_avg_shift');
  
 % vol_avg_noshift_down = vol_avg_noshift(ind_X_down,ind_Y_down,ind_Z_down,:);
%  vol_avg_shift_AVG_down = vol_avg_shift_AVG(ind_X_down,ind_Y_down,ind_Z_down,:);
%  vol_avg_noshift_AVG_down = vol_avg_noshift_AVG(ind_X_down,ind_Y_down,ind_Z_down,:);
  
  
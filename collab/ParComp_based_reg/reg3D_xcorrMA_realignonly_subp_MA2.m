
%reg3D_xcorrMA.m
%%now do 3D corregistration
function [vol_avg_shift, Mask_postreg] = reg3D_xcorrMA_realignonly_subp( vol_avg_USE2, ShiftMat, Nupsamp_XYZ);
%returns shifted volume stack, and single volume mask of pixels with
%successful interpolation in every volume of the stack

if nargin < 3
    Nupsamp_XYZ = [1 1 1];
end


%Nfromedge should be all even numbers


%WIN_IM = [1 128; 1 size(vol_avg,2); 1 (size(vol_avg,3)-1)];

Nvol = size(vol_avg_USE2,4);


Lx2 = size(vol_avg_USE2,1);
Ly2 = size(vol_avg_USE2,2);
Lz2 = size(vol_avg_USE2,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now realign using estimated shifts

  display('realigning volumes'); tic

  
  dS_X = 1/Nupsamp_XYZ(1);
  dS_Y = 1/Nupsamp_XYZ(2);
  dS_Z = 1/Nupsamp_XYZ(3);
  
  %X0 = 1:Lx2;
  %Y0 = 1:Ly2;
  %Z0 = 1:Lz2;
  
  %X1a = 1:dS_X:Lx2;
  %Y1a = 1:dS_Y:Ly2;
  %Z1a = 1:dS_Z:Lz2;
  
  Y0 = 1:Lx2;
  X0 = 1:Ly2;
  Z0 = 1:Lz2;
  
  Y1a = 1:dS_X:Lx2;
  X1a = 1:dS_Y:Ly2;
  Z1a = 1:dS_Z:Lz2;
  
  Lx3 = length(X1a);
  Ly3 = length(Y1a);
  Lz3 = length(Z1a);
  
  
[MeshX0, MeshY0, MeshZ0] = meshgrid(X0, Y0, Z0);
%[MeshX0, MeshY0, MeshZ0] = meshgrid(Y0, X0, Z0);

  
  %for subpixel version, do interpolation using interp3:
  METHOD = 'linear';
%vol_avg_shift = zeros(Lx3,Ly3,Lz3,Nvol);
vol_avg_shift = zeros(Ly3,Lx3,Lz3,Nvol);
for count = 1:Nvol
    fprintf([num2str(count),'\n']);
%    X1 = X1a + ShiftMat(count,1);
%    Y1 = Y1a + ShiftMat(count,2);
%    Z1 = Z1a + ShiftMat(count,3);

    X1 = X1a + ShiftMat(count,1); %NEW MA May 2012 reverse col 2 and col 1
    Y1 = Y1a + ShiftMat(count,2);
    Z1 = Z1a + ShiftMat(count,3);
    NEWFIX_MA = 2;
    
    OLD = 0;
    if OLD == 1
        %cludge for forcing interpolation at borders: currently, if X0 = 15, X1 = 15.000,
        %then it won't do it...
        X1(end) = X1(end) - .0001;
        Y1(end) = Y1(end) - .0001;
        Z1(end) = Z1(end) - .0001;
        X1(1) = X1(1) + .0001;
        Y1(1) = Y1(1) + .0001;
        Z1(1) = Z1(1) + .0001;
    end
    
    
    [MeshX1, MeshY1, MeshZ1] = meshgrid(X1, Y1, Z1);

%    [MeshX1, MeshY1, MeshZ1] = meshgrid(Y1, X1, Z1);
%    [MeshY1, MeshX1, MeshZ1] = meshgrid(Y1, X1, Z1);

%size(Y1)
%size(X1)
%    size(squeeze(vol_avg_USE2(:,:,:,count)))
%size(MeshX0)
%size(MeshX1)
    vol_avg_shift(:,:,:,count) = interp3(MeshX0,MeshY0,MeshZ0,squeeze(vol_avg_USE2(:,:,:,count)),MeshX1,MeshY1,MeshZ1, METHOD, -1);
end

size(vol_avg_shift)
size(vol_avg_USE2)

display('done realigning volumes'); toc

%%% now send back a mask volume of good pixels: 
Mask_postreg = squeeze(min(vol_avg_shift,[],4) > -1);





  
  
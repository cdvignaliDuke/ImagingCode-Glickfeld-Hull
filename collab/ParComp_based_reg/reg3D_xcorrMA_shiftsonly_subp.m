
%reg3D_xcorrMA.m
%%now do 3D corregistration
function [ShiftMat, ShiftMat000] = reg3D_xcorrMA_shiftsonly_subp( vol_avg, vol_avg_mean, WIN_IM, Nupsamp, Nfromedge, Nvolbin);
%interpolate the cost function to estimate subpixel shifts.
%Nfromedge should be all even numbers


%WIN_IM = [1 128; 1 size(vol_avg,2); 1 (size(vol_avg,3)-1)];
%Nupsamp = 4;
%Nfromedge = [4 4 3]; 



if nargin < 6
    Nvolbin = 1;
end
if nargin < 5
    Nfromedge = [4 4 3];
end
if nargin < 4
    Nupsamp = 4;
end
if nargin < 3
    WIN_IM = [1 30; 1 30; 1 8];
end

Nvol0 = size(vol_avg,4);
if Nvolbin == 1
%    vol_avg = vol_avg0;
%    clear('vol_avg0');
    Nvol = Nvol0;
else
    Nvol = floor(Nvol0./Nvolbin);
end


Lx = size(vol_avg,1);
Ly = size(vol_avg,2);
Lz = size(vol_avg,3);


NfromedgeX = Nfromedge(1);
NfromedgeY = Nfromedge(2);
NfromedgeZ = Nfromedge(3); %starting pixel in Z


ind_X_Im = WIN_IM(1,1):WIN_IM(1,2);
ind_Y_Im = WIN_IM(2,1):WIN_IM(2,2);
ind_Z_Im = WIN_IM(3,1):WIN_IM(3,2);

ind_X = (WIN_IM(1,1)-1+NfromedgeX):(WIN_IM(1,2)+1-NfromedgeX);
ind_Y = (WIN_IM(2,1)-1+NfromedgeY):(WIN_IM(2,2)+1-NfromedgeY);
ind_Z = (WIN_IM(3,1)-1+NfromedgeZ):(WIN_IM(3,2)+1-NfromedgeZ);

if Nupsamp == 1
  %%%version of 3d motion algorithm without subpixel
  Ref_Im = vol_avg_mean(ind_X,ind_Y,ind_Z);

  display('calculating 3D xcorr, and shift values'); tic
  display('vol #:');
  ShiftMat0 = ones(Nvol0,1)*Nfromedge;
  for count= 1:Nvol
      if Nvolbin == 1
          fprintf([num2str(count),'\n']);
          Im = vol_avg(ind_X_Im,ind_Y_Im,ind_Z_Im,count);
          C = normxcorr3(Ref_Im,Im,'valid');
          ijk_max = findn(C == max(max(max(C))));
          ShiftMat0(count,:) = ijk_max;
      elseif Nvolbin > 1
          fprintf([num2str((count-1)*Nvolbin+1),'\n']);
          Im = squeeze(mean(vol_avg(ind_X_Im,ind_Y_Im,ind_Z_Im,(count-1)*Nvolbin + [1:Nvolbin]),4));
          C = normxcorr3(Ref_Im,Im,'valid');
          ijk_max = findn(C == max(max(max(C))));
          ShiftMat0((count-1)*Nvolbin + [1:Nvolbin],:) = ones(Nvolbin,1)*ijk_max;
      end
  end
  
  display('done calculating 3D xcorr, and shift values'); toc
 % vol_avg_USE2 = vol_avg;
 % clear('vol_avg');
  
elseif Nupsamp > 1
  
  
  %%now do version of everhing with subpixel, use same upsampling in
  %all axes
  if Nupsamp == 1
    Ntimes = 0;
  elseif Nupsamp == 2
    Ntimes = 1;
  elseif Nupsamp == 4
    Ntimes = 2;
  end
%  METHOD = 'linear';
  METHOD = 'cubic';
  

X0vec = [(-NfromedgeX+1):(NfromedgeX-1)];
Y0vec = [(-NfromedgeY+1):(NfromedgeY-1)];
Z0vec = [(-NfromedgeZ+1):(NfromedgeZ-1)];

[MeshX0, MeshY0, MeshZ0] = meshgrid(X0vec, Y0vec, Z0vec);


dS2 = 1./Nupsamp;

Xvec = [(-NfromedgeX+1):dS2:(NfromedgeX-1)];
Yvec = [(-NfromedgeY+1):dS2:(NfromedgeY-1)];
Zvec = [(-NfromedgeZ+1):dS2:(NfromedgeZ-1)];
[MeshX, MeshY, MeshZ] = meshgrid(Xvec, Yvec, Zvec);

 Ref_Im = vol_avg_mean(ind_X,ind_Y,ind_Z);

  display('calculating 3D xcorr, and shift values'); tic
  display('vol #:');
  ShiftMat0 = ones(Nvol0,1)*Nfromedge;
  ShiftMat00 = ones(Nvol0,1)*Nfromedge;
  for count= 1:Nvol
      if Nvolbin == 1
          fprintf([num2str(count),'\n']);
          Im = vol_avg(ind_X_Im,ind_Y_Im,ind_Z_Im,count);
          C = normxcorr3(Ref_Im,Im,'valid');
          ijk_max = findn(C == max(max(max(C))));
          ijk_max = ijk_max(1,:);
          
         
          ShiftMat0(count,:) = ijk_max;
            %now interpolate this: 
          C2 = interp3(MeshX0,MeshY0,MeshZ0,C,MeshX,MeshY,MeshZ,METHOD);
          ijk_max2 = findn(C2 == max(max(max(C2))));;
          ShiftMat00((count-1)*Nvolbin + [1:Nvolbin],:) = ones(Nvolbin,1)*[MeshX(ijk_max2(1),ijk_max2(2),ijk_max2(3)), ...
              MeshY(ijk_max2(1),ijk_max2(2),ijk_max2(3)), ...
              MeshZ(ijk_max2(1),ijk_max2(2),ijk_max2(3))];

      elseif Nvolbin > 1
          fprintf([num2str((count-1)*Nvolbin+1),'\n']);
          Im = squeeze(mean(vol_avg(ind_X_Im,ind_Y_Im,ind_Z_Im,(count-1)*Nvolbin + [1:Nvolbin]),4));
          C = normxcorr3(Ref_Im,Im,'valid');
          %hi = 1
          %pause(10)
          ijk_max = findn(C == max(max(max(C))));
          ShiftMat0((count-1)*Nvolbin + [1:Nvolbin],:) = ones(Nvolbin,1)*ijk_max;
            %now interpolate this: 
          C2 = interp3(MeshX0,MeshY0,MeshZ0,C,MeshX,MeshY,MeshZ,METHOD);
          ijk_max2 = findn(C2 == max(max(max(C2))));;
          ShiftMat00((count-1)*Nvolbin + [1:Nvolbin],:) = ones(Nvolbin,1)*[MeshX(ijk_max2(1),ijk_max2(2),ijk_max2(3)), ...
              MeshY(ijk_max2(1),ijk_max2(2),ijk_max2(3)), ...
              MeshZ(ijk_max2(1),ijk_max2(2),ijk_max2(3))];
              
      end
  end
  
  display('done calculating 3D xcorr, and shift values'); toc
 % vol_avg_USE2 = vol_avg;
 % clear('vol_avg');
 
 % vol_avg_USE2 = vol_avg_interp;
 % clear('vol_avg');
  end

%%now find the shift from baseline: 
ShiftMat = zeros(size(ShiftMat0));
ShiftMat_Cent = Nfromedge;
%ShiftMat(:,1) = ShiftMat0(:,1) - median(ShiftMat0(:,1));
%ShiftMat(:,2) = ShiftMat0(:,2) - median(ShiftMat0(:,2));
%ShiftMat(:,3) = ShiftMat0(:,3) - median(ShiftMat0(:,3));

ShiftMat(:,1) = -1*(ShiftMat0(:,1) - ShiftMat_Cent(1)); 
ShiftMat(:,2) = -1*(ShiftMat0(:,2) - ShiftMat_Cent(2));
ShiftMat(:,3) = -1*(ShiftMat0(:,3) - ShiftMat_Cent(3));

if Nupsamp == 1
    ShiftMat000 = ShiftMat;
else
    ShiftMat000 = ShiftMat00;
end

%ShiftMat_noupsamp = ShiftMat./Nupsamp;



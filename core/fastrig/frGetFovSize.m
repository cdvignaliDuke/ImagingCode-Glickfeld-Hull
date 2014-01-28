function fov_siz = frGetFovSize(zoom_v,obj_mag);
%FRGETFOVSIZE
% FOV_SIZ = FRGETFOVSIZE(ZOOM_V,OBJ_MAG);

slope = 205.754; % um / volt

fov_siz = round(slope*zoom_v * 16 /obj_mag/10)*10;

return;

% calibration code. done by Mark H using 16x objective.

umPerL = [ repmat(1000/25, [1 3]) repmat(1000/105, [1 8])];
tZoom = [   4 3 2 ...
            2   1   0.75 0.50 0.25  0.1 0.05    0.02];
nLines = [  8 8 7 ...
            10  10  10   10   2     2   1       1  ];

nPix =  [   97.9 139.0 171.2 ...
            59.3 114.8 152.05 225.9 ... 
            90.3 204.5 156.8  203.1]


umPerPix = (nLines .*umPerL) ./ nPix;

b = regress(umPerPix(:)*256,tZoom(:));

clf
plot(tZoom(:),tZoom(:)*b,'-k')
hold on;
plot(tZoom(:),umPerPix(:)*256,'dk','markerfacecolor','w');

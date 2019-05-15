FSAV_attentionV1
exptBx = expt;
clear expt
FSAV_V1_100ms_naive
exptNv = expt;

%% Behavior
nExp = size(exptBx,2);
FOV_xum_bx = nan(1,nExp);
FOV_yum_bx = nan(1,nExp);
z_bx = nan(1,nExp);
for iexp = 1:nExp
    d = str2double(exptBx(iexp).date);
    obj = exptBx(iexp).obj;
    zoom = exptBx(iexp).zoom;
    [x_um,y_um] = scalebarCalib(d,obj);
    FOV_xum_bx(iexp) = x_um./zoom;
    FOV_yum_bx(iexp) = y_um./zoom;
    z_bx(iexp) = exptBx(iexp).z;
end

%% Naive
nExp = size(exptNv,2);
FOV_xum_nv = nan(1,nExp);
FOV_yum_nv = nan(1,nExp);
z_nv = nan(1,nExp);
for iexp = 1:nExp
    d = str2double(exptNv(iexp).date);
    obj = exptNv(iexp).obj;
    zoom = exptNv(iexp).zoom;
    [x_um,y_um] = scalebarCalib(d,obj);
    FOV_xum_nv(iexp) = x_um./zoom;
    FOV_yum_nv(iexp) = y_um./zoom;
    z_nv(iexp) = exptNv(iexp).z;
end

%% Summary

FOV_xum = mean(cat(2,FOV_xum_bx,FOV_xum_nv));
FOV_xumErr = ste(cat(2,FOV_xum_bx,FOV_xum_nv),2);
FOV_yum = mean(cat(2,FOV_yum_bx,FOV_yum_nv));
FOV_yumErr = ste(cat(2,FOV_yum_bx,FOV_yum_nv),2);
fprintf('X: %s +/- %s, Y: %s +/- %s\n', num2str(FOV_xum),num2str(FOV_xumErr),...
    num2str(FOV_yum),num2str(FOV_yumErr))
fprintf('X range: %s - %s\n',num2str(min(cat(2,FOV_xum_bx,FOV_xum_nv))),...
    num2str(max(cat(2,FOV_xum_bx,FOV_xum_nv))))
fprintf('Y range: %s - %s\n',num2str(min(cat(2,FOV_yum_bx,FOV_yum_nv))),...
    num2str(max(cat(2,FOV_yum_bx,FOV_yum_nv))))

z = mean(cat(2,z_bx,z_nv));
zErr = ste(cat(2,z_bx,z_nv),2);
zRange = [min(cat(2,z_bx,z_nv)),max(cat(2,z_bx,z_nv))];
fprintf('Z: %s +/- %s, range:%s\n',num2str(z), num2str(zErr), num2str(zRange))
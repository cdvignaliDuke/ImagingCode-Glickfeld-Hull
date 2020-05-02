function [reg2ref, reg2ref_dfof] = splitTransform(ref, reg, data_dfof_max_reg)
%% Sectioning
    refTL = ref(1:512/2-1, 1:796/2-1);
    refTR = ref(1:512/2-1, 796/2:end);
    refBL = ref(512/2:end, 1:796/2-1);
    refBR = ref(512/2:end, 796/2:end);
    regTL = reg(1:512/2-1, 1:796/2-1);
    regTR = reg(1:512/2-1, 796/2:end);
    regBL = reg(512/2:end, 1:796/2-1);
    regBR = reg(512/2:end, 796/2:end);
    data_dfof_max_regTL = data_dfof_max_reg(1:512/2-1, 1:796/2-1);
    data_dfof_max_regTR = data_dfof_max_reg(1:512/2-1, 796/2:end);
    data_dfof_max_regBL = data_dfof_max_reg(512/2:end, 1:796/2-1);
    data_dfof_max_regBR = data_dfof_max_reg(512/2:end, 796/2:end);
    
%% Transforming
% Top Left
    sz_target  = size(regTL);
    [input_points, base_points] = cpselect(double(regTL),double(refTL),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    regTL2refTL = imtransform(double(regTL),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    regTL2refTL_dfof = imtransform(double(data_dfof_max_regTL),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

% Top Right
    sz_target  = size(regTR);
    [input_points, base_points] = cpselect(double(regTR),double(refTR),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    regTR2refTR = imtransform(double(regTR),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    regTR2refTR_dfof = imtransform(double(data_dfof_max_regTR),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    
% Bottom Left
    sz_target  = size(regBL);
    [input_points, base_points] = cpselect(double(regBL),double(refBL),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    regBL2refBL = imtransform(double(regBL),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    regBL2refBL_dfof = imtransform(double(data_dfof_max_regBL),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

% Bottom Right
    sz_target  = size(regBR);
    [input_points, base_points] = cpselect(double(regBR),double(refBR),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    regBR2refBR = imtransform(double(regBR),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    regBR2refBR_dfof = imtransform(double(data_dfof_max_regBR),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

%% Combining
    reg2ref = [regTL2refTL, regTR2refTR; regBL2refBL, regBR2refBR];
    reg2ref_dfof = [regTL2refTL_dfof, regTR2refTR_dfof; regBL2refBL_dfof, regBR2refBR_dfof];   
    
    
    
    
% Original    
%     sz_target  = size(reg);
%     [input_points, base_points] = cpselect(double(reg),double(ref),'Wait', true);
%     mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
%     reg2ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
%     reg2ref_dfof = imtransform(double(data_dfof_max_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
%     
    
    
end

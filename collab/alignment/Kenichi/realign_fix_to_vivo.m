transform_fix_to_vivo = vivo_cells/fixed_cells;
transform_vivo_to_fix = fixed_cells/vivo_cells;

dim=size(vivo);
dim2=size(fixed);
fixed_aligned=zeros(dim(1),dim(2),dim(3));

for x=1:dim(1)
    x
    for y=1:dim(2)
        for z=1:dim(3)
            XYZ=round(transform_vivo_to_fix * [x;y;z;1]);
            if XYZ(1) < 1 | XYZ(2) < 1 | XYZ(3) < 1 | XYZ(1) > dim2(1) | XYZ(2) > dim2(2) | XYZ(3) > dim2(3)
                fixed_aligned(x,y,z)=0;
            else
                fixed_aligned(x,y,z)=fixed(XYZ(1),XYZ(2),XYZ(3));
            end
        end
    end
end

clear XYZ ans dim dim2 x y z
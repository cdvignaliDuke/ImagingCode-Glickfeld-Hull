function angle = write_angle_map(th, fname)

dim=size(th);
temp=ones(dim(1),dim(2));

angle=hsv2rgbKO(th, temp,temp);

imwrite(angle, [fname, 'angle.tif']);

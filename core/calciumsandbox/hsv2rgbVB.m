function rgb = hsv2rgbVB (hue, saturation, lightness) 
% DEPRECATED

if hue < 0
    hue = 0;
else
    if hue > 1;
        hue = 1;
    end
end

if saturation < 0
    saturation = 0;
else
    if saturation > 1;
        saturation = 1;
    end
end

if lightness < 0
    lightness = 0;
else
    if lightness > 1;
        lightness = 1;
    end
end

%%%%% TODO

v = ones(size(hue));

rgb = zeros([size(hue),3]);
rgb = cat(v,hue,0*v);

sel = find(hue<0.25);
red = ones(

if hue < 0.25
    red = (1-hue*4)+(hue*4);
else
    
end


if hue < 0.25
    rgb = red * (1-hue*4) + yellow * hue*4;
else
    if hue < 0.5
        hue = hue - 0.25;
        rgb = yellow * (1-hue*4) + green * hue*4;
    else if hue < 0.75
            hue = hue -0.5;
            rgb = green * (1-hue*4) + blue * hue*4;
        else
            hue = hue -0.75;
            rgb = blue * (1-hue*4) + red * hue*4;
        end
    end
end
rgb= rgb.*saturation + white.*(1-saturation);
rgb= rgb.*lightness;


                
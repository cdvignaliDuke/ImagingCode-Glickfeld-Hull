function in = imAdaptHistEqJavaRevBorders(in, win_size)
%IMADAPTHISTEQJAVAREVBORDERS
% IN = IMADAPTHISTEQJAVAREVBORDERS(IN, WIN_SIZE)
%
% See STACKLOCALCONTRASTADAPTATION

    dim=size(in);
    in=[fliplr(in),in,fliplr(in)];
    in=[flipud(in);in;flipud(in)];
    in = imAdaptHistEqJava(in, win_size);
    in=in(dim(1)+1:dim(1)*2,dim(2)+1:dim(2)*2);
end
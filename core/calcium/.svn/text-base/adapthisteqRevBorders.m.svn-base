
function in = adapthisteqRevBorders(in, Ntiles_x, Ntiles_y, cliplimit)
    if nargin < 3
        cliplimit=0.02;
    end
    dim=size(in);
    in=[fliplr(in),in,fliplr(in)];
    in=[flipud(in);in;flipud(in)];
    in = adapthisteq(in,'clipLimit',cliplimit,'numTiles',[Ntiles_x Ntiles_y]); % MATLAB image processing toolbox
    in=in(dim(1)+1:dim(1)*2,dim(2)+1:dim(2)*2);
end
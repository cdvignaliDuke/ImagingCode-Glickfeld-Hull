function [frame]=resampleframe(raw,xii,width)

[nrows,ncols]=size(raw);

frame = zeros(nrows,width,'uint16');

for icol = 1:width
    sel = find(xii==icol);
    if length(sel)>1
        frame(:,icol) = 2048-sum(raw(:,sel),2)/sum(sel);
%        frame(:,icol) = 2048-raw(:,sel(1));
    else
        frame(:,icol) = 2048-raw(:,sel(1));
    end
end

return;

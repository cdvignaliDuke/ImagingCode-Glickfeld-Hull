function map=ezCellMap(index, labelimg, cell_to_exclude, show_flag)

% Simple code to generate a map of "index"
% index: Nx1 or 1xN, indices of N cells
%
%   05.12.2008.     Kenichi Ohki

if nargin<3
    cell_to_exclude=[];
end

if nargin<4
    show_flag=1;
end

[Nx,Ny]=size(labelimg);
Ncells=max(labelimg(:));
map=zeros(Nx,Ny);
for i=1:Ncells
    if isempty(find(cell_to_exclude==i))
        map(find(labelimg==i))=index(i);
    end
end
if show_flag==1
    figure;
    imagesc(map);
    colorbar;
end

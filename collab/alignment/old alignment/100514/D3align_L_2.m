function [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(cells1, cells2, coordinate_pairs,list_pairs_in,fig_flag,tsh_dist,num_of_manual_pairs)
%[transform_coeff transform_coeff list_pairs_out]=D3align_L_2('invivo_rotated_um.txt', 'invitro_um.txt','invivo_invitro_pairs.txt');

% cells1 -  three columm matrix, containing coordinates of first cell set OR
%           file name of ascii file containing matrix
% cells2 -  three columm matrix, containing coordinates of second cell set
%           OR file name of ascii file containing matrix
% coordinate_pairs - coordinates of corresponding points in first and
%                   second sets. Text file or matrix. Format: odd lines - coordinates from
%                   first set, even line - corresponding coordinates from second set
% list_pairs_in - indexes of cell pairs. Text file or matrix. Format -
%                   first column - cell index from first set, second column - index of
%                   corresponding cell from second set.

% one of coordinate_pairs or list_pairs_in has to be specified
if nargin<6 || isempty(tsh_dist)
    tsh_dist=8; %max distance between paired cells (in pixels)
end
if nargin<7 || isempty(num_of_manual_pairs)
    num_of_manual_pairs=5; %number of pairs required to start auto alignment
end
%num_of_manual_pairs=4;
%tsh_dist=7;

transform_coeff(1:12)=0;
list_pairs_out=[];
dist_list=[];
inverse_transform_coeff(1:12)=0;

if nargin<5 || isempty(fig_flag)
    fig_flag=1;
end

%read file 1
if ischar(cells1) % if filename
    p_1=load(cells1,'-ascii');
else
    p_1=cells1;
end


%read file 2
if ischar(cells2) % if filename
    p_2=load(cells2,'-ascii');
else
    p_2=cells2;
end
p_1=p_1(:,1:3);
p_2=p_2(:,1:3);


% create corresponding coordinates pairs
if nargin<3 % no initial pairs
    display('Requires initial pair list');

    transform_coeff(1)=mean(p_1(:,1))-mean(p_2(:,1));
    transform_coeff(5)=mean(p_1(:,2))-mean(p_2(:,2));
    transform_coeff(9)=mean(p_1(:,3))-mean(p_2(:,3));
    transform_coeff(2)=1;
    transform_coeff(7)=1;
    transform_coeff(12)=1;

    h=drawCells(p_1,p_2,transform_coeff);
    
%drawCells(p_1,p_2,transform_coeff);
    %input pairs
    list_pairs=[];
    answer=requestPair;
    if length(answer)>1
        list_pairs=[list_pairs; [str2num(answer{1}) str2num(answer{2})]];
        close(h);
        list_pairs
        [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],list_pairs,fig_flag,tsh_dist,num_of_manual_pairs);
    end
    return
end


% some initial pairs available
if ~isempty(coordinate_pairs) %coordiante pairs available
    if ischar(coordinate_pairs)
        p_p=load(coordinate_pairs,'-ascii');
    else
        p_p=coordinate_pairs;
    end
    p_p=p_p(:,1:3);
    p_p_1=p_p(1:2:end,:);
    p_p_2=p_p(2:2:end,:);

    if size(p_p_1,1)<num_of_manual_pairs
        display('Requires initial pair list with more than 4 pairs');
        if size(p_p_1,1)<3 % one or two pairs - shift to mean of pairs
            transform_coeff(1)=mean(p_p_1(:,1))-mean(p_p_2(:,1));
            transform_coeff(5)=mean(p_p_1(:,2))-mean(p_p_2(:,2));
            transform_coeff(9)=mean(p_p_1(:,3))-mean(p_p_2(:,3));
            transform_coeff(2)=1;
            transform_coeff(7)=1;
            transform_coeff(12)=1;
        else
            % for more than two pairs - approximate simple transformation using regression
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,1)];
            transform_coeff(1:2)=regress(p_p_1(:,1),Xr);
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,2)];
            transform_coeff([5 7])=regress(p_p_1(:,2),Xr);
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,3)];
            transform_coeff([9 12])=regress(p_p_1(:,3),Xr);
        end

        h=drawCells(p_1,p_2,transform_coeff);

        %input pairs
        answer=requestPair;
        new_coordinate_pairs=[p_p; [p_1(str2num(answer{1}),1) p_1(str2num(answer{1}),2) p_1(str2num(answer{1}),3)]];
        new_coordinate_pairs=[new_coordinate_pairs; [p_2(str2num(answer{2}),1) p_2(str2num(answer{2}),2) p_2(str2num(answer{2}),3)]];

        close(h);
        new_coordinate_pairs
        [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, new_coordinate_pairs,fig_flag,tsh_dist,num_of_manual_pairs);

        return
    end

elseif nargin>3 % list pairs avalable
    if ischar(list_pairs_in)
        p_p_ind=load(list_pairs_in,'-ascii');
    else
        p_p_ind=list_pairs_in;
    end
    p_p_1=p_1(p_p_ind(:,1),:);
    p_p_2=p_2(p_p_ind(:,2),:);

    if size(p_p_ind,1)<num_of_manual_pairs
        display('Requires initial pair list with more than 4 pairs');
        if size(p_p_ind,1)<3 % one or two pairs - shift to mean of pairs
            transform_coeff(1)=mean(p_p_1(:,1))-mean(p_p_2(:,1));
            transform_coeff(5)=mean(p_p_1(:,2))-mean(p_p_2(:,2));
            transform_coeff(9)=mean(p_p_1(:,3))-mean(p_p_2(:,3));
            transform_coeff(2)=1;
            transform_coeff(7)=1;
            transform_coeff(12)=1;

        else
            % for more than two pairs - approximate simple transformation using regression
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,1)];
            transform_coeff(1:2)=regress(p_p_1(:,1),Xr);
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,2)];
            transform_coeff([5 7])=regress(p_p_1(:,2),Xr);
            Xr = [ones(size(p_p_2,1),1) p_p_2(:,3)];
            transform_coeff([9 12])=regress(p_p_1(:,3),Xr);
        end
        h=drawCells(p_1,p_2,transform_coeff);

        %input pairs
        answer=requestPair;
        list_pairs=[list_pairs_in; [str2num(answer{1}) str2num(answer{2})]];
        close(h);
        list_pairs
        [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],list_pairs,fig_flag,tsh_dist,num_of_manual_pairs);

        return
    end
end



% approximate transformation using regression
Xr = [ones(size(p_p_2,1),1) p_p_2(:,1) p_p_2(:,2) p_p_2(:,3)];
warning off all
transform_coeff(1:4)=regress(p_p_1(:,1),Xr);
transform_coeff(5:8)=regress(p_p_1(:,2),Xr);
transform_coeff(9:12)=regress(p_p_1(:,3),Xr);
warning on all

% %check transformation for "perpendicularity" and stop if "not
% perpendicular"
[a0 b0 c0]=transform(transform_coeff,0,0,0);
[a1 b1 c1]=transform(transform_coeff,0,0,1);
dtr=sqrt((a0-a1).^2 +(b0-b1).^2+(b0-b1).^2);
if dtr==0
    return
end
if abs((c0-c1)/dtr)>0.7
    return
end


% apply first transformation
[outx outy outz]=transform(transform_coeff,p_2(:,1),p_2(:,2),p_2(:,3));


% find all good pairs

sz1=size(p_1,1);
sz2=length(outx);
dx=repmat(p_1(:,1),[1 sz2])-repmat(outx',[sz1 1]);
dy=repmat(p_1(:,2),[1 sz2])-repmat(outy',[sz1 1]);
dz=repmat(p_1(:,3),[1 sz2])-repmat(outz',[sz1 1]);
d3D=sqrt((dx).^2+(dy).^2+(dz).^2); 



list_pairs_out=[];
dist_list=[];
i_list=[];
for p=1:sz1
    [C,I] = min(d3D(:));
    if C(1)<tsh_dist
        [i,j]=ind2sub([sz1 sz2],I(1));
        list_pairs_out=[list_pairs_out; [i j]];
        d3D(i,:)=inf;
        d3D(:,j)=inf;
        dist_list=[dist_list C(1)];
        i_list=[i_list i];
    else
        break
    end
end
if fig_flag
display (['Number of pairs = ' num2str(size(list_pairs_out,1)) ', mean distance = ' num2str(mean(dist_list))]);
end

% if find extra pairs then search for more
if size(list_pairs_out,1)>size(p_p_1,1)
    [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],list_pairs_out,fig_flag,tsh_dist,num_of_manual_pairs);
else
    
% calculate inverse transformation here
inverse_transform_coeff=calculateRevTransform(transform_coeff);

if fig_flag
%draw cells
    drawCells(p_1,p_2,transform_coeff);
    figure(1002);
    hist(dist_list,0:tsh_dist);
end
end


function answer=requestPair()
prompt = {'Enter fist cell in pair (c):','Enter second cell in pair (m):'};
dlg_title = 'Input cells pairs';
num_lines = 1;
def = {'',''};
options.WindowStyle='normal';
answer = inputdlg(prompt,dlg_title,num_lines,def,options);

function h=drawCells(p1,p2,tr_coeffs)

%h=figure;
h=figure(1000);
scatter3(p1(:,1),p1(:,2),p1(:,3),'c','filled');
for n=1:length(p1(:,1))
    text(p1(n,1)+4,p1(n,2)+4,p1(n,3)+4,num2str(n),'Color','c');
end
hold on

[outx outy outz]=transform(tr_coeffs,p2(:,1),p2(:,2),p2(:,3));

scatter3(outx,outy,outz,'m','filled');
for n=1:length(p2(:,1))
    text(outx(n)-3,outy(n)-3,outz(n)-3,num2str(n),'Color','m');
end
hold off;


function pairs=findPotentialPairs(p_1, p_2)

p_1=p_1(:,1:3);
p_2=p_2(:,1:3);


le1=size(p_1,1);
le2=size(p_2,1);

inter_dist_1=zeros(le1);
inter_dist_2=zeros(le2);

z_comp_1=zeros(le1);
z_comp_2=zeros(le2);

xy_comp_1=zeros(le1);
xy_comp_2=zeros(le2);


% calculate pairwise distance, z component and xy component for set1
for i=1:le1
    dist=sqrt((p_1(:,1)-p_1(i,1)).^2 + (p_1(:,2)-p_1(i,2)).^2 + (p_1(:,3)-p_1(i,3)).^2);
    inter_dist_1(:,i)=dist;

    dist=abs(p_1(:,3)-p_1(i,3));
    z_comp_1(:,i)=dist;

    dist=sqrt((p_1(:,1)-p_1(i,1)).^2 + (p_1(:,2)-p_1(i,2)).^2);
    xy_comp_1(:,i)=dist;

end

% calculate pairwise distance, z component and xy component for set2
for i=1:le2
    dist=sqrt((p_2(:,1)-p_2(i,1)).^2 + (p_2(:,2)-p_2(i,2)).^2 + (p_2(:,3)-p_2(i,3)).^2);
    inter_dist_2(:,i)=dist;

    dist=abs(p_2(:,3)-p_2(i,3));
    z_comp_2(:,i)=dist;

    dist=sqrt((p_2(:,1)-p_2(i,1)).^2 + (p_2(:,2)-p_2(i,2)).^2);
    xy_comp_2(:,i)=dist;
end


%{
W=zeros(le1,le2);
[max_d_2, max_i_2]=sort(inter_dist_2(:));

for t=1:length(max_d_2)
    [fx1, fy1]=find(inter_dist_2==max_d_2(t));
    [fx2, fy2]=find(inter_dist_1<(max_d_2(t)*1.5+5) & inter_dist_1>(max_d_2(t)/1.5-5));
    %      W(unique(fx2),unique(fx1))=W(unique(fx2),unique(fx1))+1;
    for t1=1:length(fx1)
        for t2=1:length(fx2)
            W(fx2(t2),fx1(t1))=W(fx2(t2),fx1(t1))+1;
        end
    end
end


WZ=zeros(le1,le2);
[max_d_2, max_i_2]=sort(z_comp_2(:));

for t=1:length(max_d_2)
    [fx1, fy1]=find(z_comp_2==max_d_2(t));
    [fx2, fy2]=find(xy_comp_1<(max_d_2(t)*1.5+5) & xy_comp_1>(max_d_2(t)/1.5-5));
    %      W(unique(fx2),unique(fx1))=W(unique(fx2),unique(fx1))+1;
    for t1=1:length(fx1)
        for t2=1:length(fx2)
            WZ(fx2(t2),fx1(t1))=WZ(fx2(t2),fx1(t1))+1;
        end
    end
end


WXY=zeros(le1,le2);
[max_d_2, max_i_2]=sort(xy_comp_2(:));

for t=1:length(max_d_2)
    [fx1, fy1]=find(xy_comp_2==max_d_2(t));
    [fx2, fy2]=find(z_comp_1<(max_d_2(t)*1.5+5) & z_comp_1>(max_d_2(t)/1.5-5));
    %      W(unique(fx2),unique(fx1))=W(unique(fx2),unique(fx1))+1;
    for t1=1:length(fx1)
        for t2=1:length(fx2)
            WXY(fx2(t2),fx1(t1))=WXY(fx2(t2),fx1(t1))+1;
        end
    end
end
%}

% find cell-to-cell segmentd with similar distance, and  components
% increment correspondence matix if segments are similar
WW=zeros(le1,le2);
for t1=1:le1
    for t1a=t1+1:le1
    [fx2, fy2]=find(inter_dist_2<(inter_dist_1(t1,t1a)*1.5+6) & inter_dist_2>(inter_dist_1(t1,t1a)/1.5-6) & ...
        z_comp_2<(xy_comp_1(t1,t1a)*1.5+6) & z_comp_2>(xy_comp_1(t1,t1a)/1.5-6) & ...
        xy_comp_2<(z_comp_1(t1,t1a)*1.5+6) & xy_comp_2>(z_comp_1(t1,t1a)/1.5-6));

    for t2=1:length(fx2)
%            WW([t1 t1a],[fx2(t2) fy2(t2)])=WW([t1 t1a],[fx2(t2) fy2(t2)])+1/length(fx2);
%             WW([t1 t1a],[fx2(t2) fy2(t2)])=WW([t1 t1a],[fx2(t2) fy2(t2)])+1/max(1000,(abs(inter_dist_2(fx2(t2),fy2(t2))-inter_dist_1(t1,t1a)))+...
%             (abs(z_comp_2(fx2(t2),fy2(t2))-xy_comp_1(t1,t1a))) + ...
%             (abs(xy_comp_2(fx2(t2),fy2(t2))-z_comp_1(t1,t1a))));
        
        WW([t1 t1a],[fx2(t2) fy2(t2)])=WW([t1 t1a],[fx2(t2) fy2(t2)])+1;
    end
    end
end



% cut correspondence matrix at tsh_level
tsh_level=0.5;

ws=sort(WW(:));
tsh=ws(round(length(ws)*tsh_level));

WWW=WW*0;
WWW(WW>tsh)=1;
[p1 p2]=find(WWW);
pairs=[p1 p2];


% figure;
% imagesc(WW);
% colorbar
% 
% figure;
% imagesc(WWW);
% colorbar






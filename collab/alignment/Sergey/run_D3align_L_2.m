function [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=run_D3align_L_2(cells1, cells2, test_size,num_of_pairs,fixed_pairs)

% if pairs list is not partial, then use D3alin_L_2 directly


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

if nargin<4
    num_of_manual_pairs=5;
else
    num_of_manual_pairs=num_of_pairs;
end

p_1=p_1(:,1:3);
p_2=p_2(:,1:3);

best=0;
best_m_d=inf;
if nargin<3 || isempty(test_size)
    test_size=1000000;
end

for t=1:test_size
    rp1 = randperm(size(p_1,1));
    rp2 = randperm(size(p_2,1));

    list_pairs_in(:,1)=rp1(1:num_of_manual_pairs);
    list_pairs_in(:,2)=rp2(1:num_of_manual_pairs);
    if nargin>4
        psz=size(fixed_pairs);
        list_pairs_in(1:psz(1),:)=fixed_pairs;        
    end

    [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],list_pairs_in,0);
    warning off all
    m_d=mean(dist_list);
    warning on all

%    if m_d<3.5 && (size(list_pairs_out,1)>best || (size(list_pairs_out,1)==best && best_m_d>=m_d))
    if m_d<5 && (size(list_pairs_out,1)>best || (size(list_pairs_out,1)==best && best_m_d>=m_d))
        
        best_list_pairs_in=list_pairs_in;
        best=size(list_pairs_out,1);
        best_m_d=m_d;
        t
        D3align_L_2(p_1, p_2, [],best_list_pairs_in,1);
        drawnow;
    end
    if mod(t,500000)==0
        t
        [best best_m_d]
         drawnow;
    end
end

[transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],best_list_pairs_in,1);
function [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=run_D3align_L_2(cells1, cells2, test_size,num_of_pairs,fixed_pairs,analysis_flags)

if nargin<6 || isempty(analysis_flags)
   analysis_flags=[1 1 1]; 
end

global stopFlag;
    stopFlag=0;
    
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

if nargin<4 || isempty(num_of_pairs)
    num_of_manual_pairs=5;
else
    num_of_manual_pairs=num_of_pairs;
end

p_1=p_1(:,1:3);
p_2=p_2(:,1:3);
if analysis_flags(1)
%reduce cell-to-cell correspondence space
pairs=findPotentialPairs(p_1, p_2);
end

best=0;
best_m_d=inf;
if nargin<3 || isempty(test_size)
    test_size=1000000;
end

figure(1000);
set(1000,'KeyPressFcn',@stopA);
figure(1002);
set(1002,'KeyPressFcn',@stopA);

for t=1:test_size
    if mod(t,100)==0
        drawnow
    end
    if stopFlag
        break
    end
if analysis_flags(1)
 %renadomly permutate pair list
     rp1 =round(rand(num_of_manual_pairs,1)*(size(pairs,1)-1))+1;
     list_pairs_in(:,1)=pairs(rp1,1);
     list_pairs_in(:,2)=pairs(rp1,2);   
    
else
    % if do not use findPotentialPairs
    rp1 = randperm(size(p_1,1));
    rp2 = randperm(size(p_2,1)); 
    list_pairs_in(:,1)=rp1(1:num_of_manual_pairs);
    list_pairs_in(:,2)=rp2(1:num_of_manual_pairs);

end
    


% calculate "momentum" (not exctly) of two cell subsets
dist_x_1=p_1(list_pairs_in(:,1),1)-p_1(list_pairs_in(1,1),1);
dist_y_1=p_1(list_pairs_in(:,1),2)-p_1(list_pairs_in(1,1),2);
dist_z_1=p_1(list_pairs_in(:,1),3)-p_1(list_pairs_in(1,1),3);
m_1=sum(dist_x_1.^2+dist_y_1.^2+dist_z_1.^2);


dist_x_2=p_2(list_pairs_in(:,2),1)-p_2(list_pairs_in(1,2),1);
dist_y_2=p_2(list_pairs_in(:,2),2)-p_2(list_pairs_in(1,2),2);
dist_z_2=p_2(list_pairs_in(:,2),3)-p_2(list_pairs_in(1,2),3);
m_2=sum(dist_x_2.^2+dist_y_2.^2+dist_z_2.^2);
if analysis_flags(2)
% make sure that "momentums" are not very different
if m_1<m_2/1.2 || m_1>m_2*1.2
    continue
end
end


    
    if nargin>4 && ~isempty(fixed_pairs)
        psz=size(fixed_pairs);
        list_pairs_in(1:psz(1),:)=fixed_pairs;        
    end

    [transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],list_pairs_in,0,[],num_of_manual_pairs,analysis_flags);
    warning off all
    m_d=mean(dist_list);
    warning on all

%    if m_d<3.5 && (size(list_pairs_out,1)>best || (size(list_pairs_out,1)==best && best_m_d>=m_d))
    if m_d<5 && (size(list_pairs_out,1)>best || (size(list_pairs_out,1)==best && best_m_d>=m_d))
        
        best_list_pairs_in=list_pairs_in;
        best=size(list_pairs_out,1);
        best_m_d=m_d;
        
        D3align_L_2(p_1, p_2, [],best_list_pairs_in,1);

        
        %check transformation for "perpendicularity"
[a0 b0 c0]=transform(transform_coeff,0,0,0);
[a1 b1 c1]=transform(transform_coeff,0,0,1);
sprintf('%i %6.3f %6.3f',[t m_1/m_2 (c0-c1)/sqrt((a0-a1).^2 +(b0-b1).^2+(b0-b1).^2)])
display(' ');
        
        drawnow;
    end
    if mod(t,500000)==0
      sprintf('%i %6.3f %6.3f',  [t best best_m_d])
display(' ');      
         drawnow;
    end
end

[transform_coeff inverse_transform_coeff list_pairs_out dist_list]=D3align_L_2(p_1, p_2, [],best_list_pairs_in,1,[],num_of_manual_pairs,analysis_flags);

end

function stopA(src,evnt)
    global stopFlag;
    if evnt.Character == 'e'
        stopFlag=1;
    end
end


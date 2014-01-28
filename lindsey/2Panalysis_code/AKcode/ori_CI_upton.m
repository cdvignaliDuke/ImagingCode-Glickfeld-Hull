function OriStatKO = ori_CI_upton (OriStatKO,alpha)

Ncells = length(OriStatKO);
si=size(OriStatKO(1).norm_table);

for i=1:Ncells
    norm_table = OriStatKO(i).norm_table;
    ori_table = mean(reshape(norm_table(:,1:(si(2))),si(1),((si(2))/2),2),3);
    for j=1:si(1)
        [a,b,c] = vector_average(ori_table(j,:)./norm_table(j,si(2)));
        angles(j) = a;
    end
    OriStatKO(i).ori_vector_angle_CI = angle_confidenceAK(angles,alpha)/2;
end
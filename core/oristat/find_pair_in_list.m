function answer = find_pair_in_list (list, pair)

for i=1:size(list,1);
    if list(i,:) == pair
        answer =1;
        return;
    end
end
answer =0;

function newArray = combineIndAcrossCellArray(array)

[nexp,~] = size(array);
newArray = cell(nexp,1);
for i = 1:nexp
    newArray{i} = sum(cell2mat(array(i,:)'),1) > 0;
end

end
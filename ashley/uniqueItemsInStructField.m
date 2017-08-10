function items = uniqueItemsInStructField(dataStruct,fn)
% get a list of all unique items in a given data structure under the given
% fieldname, fn
% fn should be a string and dataStruct should only have one level

n = size(dataStruct,2);
itemsInExpt = [];
for i = 1:n
    d = eval(['dataStruct(i).' fn]);
    itemsInExpt = cat(2,itemsInExpt,d);
end
items = unique(itemsInExpt);
if iscell(items)
    if ischar(items{1})
        items = sort_nat(items);
    end
end
end
function im = maskSelect(mask,list)
% IM = MASKSELECT(MASK,LIST)

im = zeros(size(mask));

for ind = 1:length(list)
    im = im | mask==list(ind);
end

return;


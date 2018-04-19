mask = zeros(size(mask_raw));
k=1;
for i = 1:max(max(mask_raw))
    if ~isempty(mask_raw==i)
        mask(mask_raw==i) = k;
        k = k+1;
    end
end
mask_raw = mask;
mask_final = reshape(mask, 1, npw*nph);
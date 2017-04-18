function [new_mat1, new_mat2] = catIfIndexNotEmpty(cat_mat1,og_mat1,cat_mat2,og_mat2,ind,dim)  
   if all(cellfun(@isempty,ind))
       new_mat1 = cat(dim,cat_mat1,[]);
       new_mat2 = cat(dim,cat_mat2,[]);
   else
       new_mat1 = cat(dim,cat_mat1, og_mat1);
       new_mat2 = cat(dim,cat_mat2, og_mat2);
   end
end
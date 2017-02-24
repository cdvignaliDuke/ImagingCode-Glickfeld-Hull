function data_sub = subNegData(data)
%remove negative data from image stack (3D) (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);

end
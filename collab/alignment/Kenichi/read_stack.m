% Concatenate 2D images into a 3D matrix
  % Select source files.
  [FileName,PathName] = uigetfile; 

  b=zeros(256);
  DIR=dir(PathName);
  nfile=size(DIR,1)
  
  for j=[4:nfile]
    FileName=DIR(j,1).name
    a=imread([PathName FileName]);
    b=cat(3,b,a);
  end
  
  b=b(:,:,[2:size(b,3)]);
  clear DIR FileName PathName a j nfile
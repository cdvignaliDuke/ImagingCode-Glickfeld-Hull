
clear all; clc;
% Concatenate 2D images into a 3D matrix
  % Select source files.
  [FileName,PathName] = uigetfile('*.tif','Select image source'); 
  % Set a range for processing.
  prompt = {'Start:','End:'};
  dlg_title = 'Range for processing';
  answer=inputdlg(prompt,dlg_title);
  i=[str2num(answer{1}):str2num(answer{2})];
  %
  b=zeros(512);
  for j=i;                                                         % 'j' is necessary, don't replace it by 'i'
    pat = 'z(\d)(\d)(\d)';
    if (j>100), FileName=regexprep(FileName, pat, ['z' num2str(j)]);end   % regexprep (regular expression replacement)
    if (j<100)&&(j>9), FileName=regexprep(FileName, pat, ['z0' num2str(j)]);end;
    if (j<10), FileName=regexprep(FileName, pat, ['z00' num2str(j)]);end
    a=imread([PathName FileName]);
    b=cat(3,b,a);
  end
  n=size(i,2)+1;
  S=b(:,:,2:n);

% Find cells in a 3D marix.
% S is a 3D  matrix (graylevel).
% SobrdBW is a 3D  matrix (logical); '1's represent cells.
se=strel('disk',2);
Se=imerode(S,se);
Sobr=imreconstruct(Se,S);
Sobrd=imdilate(Sobr,se);
SobrdBW=imregionalmax(Sobrd);
SobrdBWclb=imclearborder(SobrdBW);

% Extract images from the 3D matrix SobrdBW.
% SobrdBW is a 3D logical matrix.
  % Select a directory for saving data.
  pause(5)
  dirName = uigetdir;
for j=[1:size(i,2)]
    QQ=SobrdBWclb(:,:,j);
    if i(j)>100, imwrite(QQ,[dirName '\z' num2str(i(j)) '_ch01.tif'],'Compression','none');end;
    if (i(j)<100)&&(i(j)>9), imwrite(QQ,[dirName '\z0' num2str(i(j)) '_ch01.tif'],'Compression','none');end;
    if i(j)<10, imwrite(QQ,[dirName '\z00' num2str(i(j)) '_ch01.tif'],'Compression','none');end;
end

%Save the 3D logical matrix
save([dirName '\BW3D'], 'SobrdBWclb')


% Find centroids of cells in a 3D logic array
% SobrdBWclb is the 3D logic array, connedted "1"s represent a cell. 
[L, num]=bwlabeln(SobrdBWclb);
s=regionprops(L,'centroid');
centroids=cat(1,s.Centroid);

% Plot centroids 
x=centroids(:,1);
y=centroids(:,2);
z=centroids(:,3);
c=[1/size(x,1):1/size(x,1):1];
figure; scatter3(x,y,z,50,c,'filled');
axis([0 512 0 512 0 512]);

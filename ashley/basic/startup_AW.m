currentdir = 'C:\Users\ashley\Documents\Repositories';
svnroot = [currentdir,'\ImagingCode-Glickfeld-Hull'];
ijroot = 'C:\Program Files\ImageJ';
addpath(genpath([currentdir,'\ImagingCode-Glickfeld-Hull']));
coreInitJavaPath(svnroot,ijroot);
coreInitMatlabPath(svnroot,ijroot);


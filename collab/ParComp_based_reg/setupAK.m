currentdir = pwd;
cd([currentdir,'\code\core\'])
svnroot = [currentdir,'\code'];
ijroot = 'D:\Program Files (x86)\ImageJ';
coreInitJavaPath(svnroot,ijroot);
coreInitMatlabPath(svnroot,ijroot);

addpath(genpath([currentdir,'\code\users\aaron\CurrentWorkingSet']));
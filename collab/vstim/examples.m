%% spatial frequency tuning
prot= readprot('prots/sftuning.prot')
makestims(prot,s,'sftuning')

%% retinotopy 
prot= readprot('prots/xpos.prot')
makestims(prot,s,'xpos');

%% preprocess pictures
sourcedir = 'C:\Documents and Settings\V1msq\My Documents\stimulation\pictures\orig\Flowers';
targetdir = 'C:\Documents and Settings\V1msq\My Documents\stimulation\pictures\proc\Flowers';
picspreprocess(sourcedir,targetdir,screen)

%% flicker pictures
[prot]= readprot('prots/pics.prot')
makestims(prot,s,'pics');
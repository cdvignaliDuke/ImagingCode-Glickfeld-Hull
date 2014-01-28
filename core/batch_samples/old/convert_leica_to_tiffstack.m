root_dir='\\zmey\storlab\data\Kenichi\mouse090509\'

% 
% data_dir='mouse090205_5'
% fname='vstim5';
% stack=readtiff([root_dir,data_dir],[],fname,1);
% writetiff(stack, [root_dir,fname]);

data_dir='mouse090509_4'
fname='vstim4';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090205_6'
fname='vstim6b';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090205_7'
fname='vstim7';
stack=readtiff([root_dir,data_dir],[],fname,1);
nframes=size(stack,3);
stackR=stack(:,:,[1:2:nframes-1]);
stackG=stack(:,:,[2:2:nframes]);
writetiff(stackR, [root_dir,fname,'_red']);
writetiff(stackG, [root_dir,fname]);


data_dir='mouse090205_8'
fname='vstim8';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090205_9'
fname='vstim9';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090205_10'
fname='vstim10';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090205_11'
fname='vstim11';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_12'
fname='vstim12';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_13'
fname='vstim13';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_14'
fname='vstim14';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_15'
fname='vstim15';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_16'
fname='vstim16';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_17'
fname='vstim17';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_18'
fname='vstim18';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_19'
fname='vstim19';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_20'
fname='vstim20';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_21'
fname='vstim21';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_22'
fname='vstim22';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_23'
fname='vstim23';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_24'
fname='vstim24';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_25'
fname='vstim25';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


data_dir='mouse090205_26'
fname='vstim26';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);


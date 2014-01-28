root_dir='E:\users\kenichi\leicatemp\'


data_dir='mouse090416_13'
fname='vstim13';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_14'
fname='vstim14';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_15'
fname='vstim15';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_16'
fname='vstim16';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_17'
fname='vstim17';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='temp'
fname='vstim18';
stack=readtiff([root_dir,data_dir],[1:1680],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_19'
fname='vstim19';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

data_dir='mouse090416_20'
fname='vstim20';
stack=readtiff([root_dir,data_dir],[],fname,1);
writetiff(stack, [root_dir,fname]);

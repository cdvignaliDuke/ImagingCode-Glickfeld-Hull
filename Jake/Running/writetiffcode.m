%SBX to Tiff
fn = 'img15_000_000';
nframes= 1000;
pn= 'D:\Jake_temp\150120img15';
pn_out = [pn 'analysis\'];

%load, baseline, register
cd = pn;
z= squeeze(sbxread(fn,0,2000));  
zz = writetiff(z,'150120img15tiff')
clear z;
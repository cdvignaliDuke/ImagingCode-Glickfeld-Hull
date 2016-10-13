clear;
day = '160314_img38';
image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
image_dest   = ['Z:\tempHoldingArea\full_size_analysis\' day]; %stores the data on crash in the lever analysis folder
bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
    mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
end
old_cd = cd; %save old cd so I can restore it later


WF_draw_ROIs_for_lever_full_size(day, image_source, image_dest);
WF_lever_Bx_DFoF_full_size;
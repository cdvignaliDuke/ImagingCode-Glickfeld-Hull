
outfname='orimap_scmr7_080322_006_p05.tif';
alpha=0.05;


% offset=45;,cmap=[[0,0,0];[0.5,0.5,0.5]];,for i=1:180, cmap=[cmap;squeeze(hsv2rgbKO(mod(180-i+1+offset,180)/180,1,1))'];,end


k=find([Oristats.p_value_ori_sel]>alpha);
k2=find([Oristats.p_value_ori_sel]<alpha | [Oristats.p_value_resp]>alpha);
map=ezCellMap(floor([Oristats.ori_vector_angle])+2, labelimg, k);
map2=ezCellMap(ones(1,length(Oristats)), labelimg, k2);
map3=map+map2;
imwrite(uint8(map3),cmap,outfname,'compression','none');

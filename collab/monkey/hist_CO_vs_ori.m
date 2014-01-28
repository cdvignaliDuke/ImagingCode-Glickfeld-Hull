function hist_CO_vs_ori (Oristats, COvalue)

alpha=0.01;
width_th=45;
width_th2=90;
OSI_th=0.3;
OSI_th2=0.15;

% statistical significance
p_resp=[Oristats.p_value_resp];
p_sel=[Oristats.p_value_sel];

% ori selectivity
OSI=[Oristats.ori_vector_tune];
ori_width=[Oristats.ori_tuning_width];
ori_width(find(isnan(ori_width)))=90;

% classify cells

SharpTuned=zeros(10,1);
BroadTuned=zeros(10,1);
Biased=zeros(10,1);
UnTuned=zeros(10,1);
NoResp=zeros(10,1);
Total=zeros(10,1);

SharpTuned2=zeros(10,1);
BroadTuned2=zeros(10,1);
Biased2=zeros(10,1);

for i=1:10
    filter=find(COvalue>=(i-1)*0.1 & COvalue <i*0.1);
    Total(i)=length(filter);
    SharpTuned(i)=length(intersect(filter,find(p_resp < alpha & p_sel < alpha & ori_width < width_th)));
    BroadTuned(i)=length(intersect(filter,find(p_resp < alpha & p_sel<alpha & ori_width >= width_th & ori_width<width_th2)));
    Biased(i)=length(intersect(filter,find(p_resp < alpha & p_sel<alpha & ori_width>=width_th2)));
    UnTuned(i)=length(intersect(filter,find(p_resp<alpha & p_sel>=alpha)));
    NoResp(i)=length(intersect(filter,find(p_resp >= alpha)));
    SharpTuned2(i)=length(intersect(filter,find(p_resp < alpha & p_sel < alpha & OSI > OSI_th)));
    BroadTuned2(i)=length(intersect(filter,find(p_resp < alpha & p_sel<alpha & OSI <= OSI_th & OSI >OSI_th2)));
    Biased2(i)=length(intersect(filter,find(p_resp < alpha & p_sel<alpha & OSI<=OSI_th2)));

end

Resp=Total-NoResp;
Tuned=Total-NoResp-UnTuned;

subplot(5,4,1);
bar(Resp./Total);
title('CO vs %Resp');

subplot(5,4,2);
bar(Tuned./Total);
title('CO vs %Tuned');

subplot(5,4,3);
bar(Tuned./Resp);
title('CO vs %Tuned of Resp');

subplot(5,4,5);
bar(SharpTuned./Total);
title('CO vs %SharpTuned (tw<45)');

subplot(5,4,6);
bar(SharpTuned./Resp);
title('CO vs %SharpTuned (tw<45) of Resp');

subplot(5,4,7);
bar(BroadTuned./Total);
title('CO vs %BroadTuned (45<tw<90)');

subplot(5,4,8);
bar(BroadTuned./Resp);
title('CO vs %BroadTuned (45<tw<90) of Resp');

subplot(5,4,9);
bar(Biased./Total);
title('CO vs %Biased (tw>90)');

subplot(5,4,10);
bar(Biased./Resp);
title('CO vs %Biased (tw>90) of Resp');

subplot(5,4,11);
bar(UnTuned./Total);
title('CO vs %UnTuned');

subplot(5,4,12);
bar(UnTuned./Resp);
title('CO vs %UnTuned of Resp');


subplot(5,4,13);
bar(SharpTuned2./Total);
title('CO vs %SharpTuned (OSI>0.3)');

subplot(5,4,14);
bar(SharpTuned2./Resp);
title('CO vs %SharpTuned (OSI>0.3) of Resp');

subplot(5,4,15);
bar(BroadTuned2./Total);
title('CO vs %BroadTuned (0.15<OSI<0.3)');

subplot(5,4,16);
bar(BroadTuned2./Resp);
title('CO vs %BroadTuned (0.15<OSI<0.3) of Resp');

subplot(5,4,17);
bar(Biased2./Total);
title('CO vs %Biased (OSI<0.15)');

subplot(5,4,18);
bar(Biased2./Resp);
title('CO vs %Biased (OSI<0.15) of Resp');

subplot(5,4,19);
bar([SharpTuned./Total,BroadTuned./Total,Biased./Total,UnTuned./Total,NoResp./Total],'stack');
title('CO vs Classification with tw');

subplot(5,4,20);
bar([SharpTuned2./Total,BroadTuned2./Total,Biased2./Total,UnTuned./Total,NoResp./Total],'stack');
title('CO vs Classification with OSI');







    

    


    













function doOriStatAK(RunInfo)

alpha = 0.05;

tc=RunInfo.Trace.tc;
nptc=RunInfo.Trace.nptc;
tcnorm=RunInfo.Trace.tcnorm;
tcnorm_np=RunInfo.Trace.tcnorm_np;
tcnorm_blank=RunInfo.Trace.tcnorm_blank;

si = size(tcnorm);
nFrames = si(1);
nreps = si(2);
ndir=si(3);
Ncells = si(4);

%Avg response tables1
onWindow = [1:length(RunInfo.fbaseWindow)]+ RunInfo.Noff;
resp_data = nanmean(tc(onWindow,:,:,:)); 
base_data = nanmean(tc(RunInfo.fbaseWindow,:,:,:));
base_data_np = nanmean(nptc(RunInfo.fbaseWindow,:,:,:));

%Avg response tables2
on_data = nanmean(tc((RunInfo.Noff+1):(RunInfo.Noff+RunInfo.Non),:,:,:)); 
off_data = nanmean(tc(1:RunInfo.Noff,:,:,:));
on_data_np = nanmean(nptc((RunInfo.Noff+1):(RunInfo.Noff+RunInfo.Non),:,:,:)); 
off_data_np = nanmean(nptc(1:RunInfo.Noff,:,:,:));

%Avg norm response table
norm_data = nanmean(tcnorm((RunInfo.Noff+1):(RunInfo.Noff+4),:,:,:)); 
norm_data_np = nanmean(tcnorm_np((RunInfo.Noff+1):(RunInfo.Noff+4),:,:,:));

data_table = squeeze(on_data);
data_table(:,(si(3)+1),:) = squeeze(nanmean(off_data,3)); 
norm_table = squeeze(norm_data);

data_table=double(data_table);
norm_table=double(norm_table);

OriStat.Ori = OriStatAK2(data_table,norm_table,alpha);

data_table_np = squeeze(on_data_np);
data_table_np(:,(si(3)+1),:) = squeeze(nanmean(off_data_np,3));
norm_table_np = squeeze(norm_data_np);
data_table_np = squeeze(nanmean(data_table_np,3));
norm_table_np = squeeze(nanmean(norm_table_np,3));
data_table_np2(:,:,1)=data_table_np;
norm_table_np2(:,:,1)=norm_table_np;
data_table_np2(:,:,2)=data_table_np;
norm_table_np2(:,:,2)=norm_table_np;

data_table_np2=double(data_table_np2);
norm_table_np2=double(norm_table_np2);

Ori_np = OriStatAK2(data_table_np2,norm_table_np2,alpha);

%plot data
tcnorm_avg_ori=squeeze(nanmean(tcnorm,2));
tcplot_tcnorm_ori=reshape(tcnorm_avg_ori(8:24,:,:),17*ndir,Ncells);

area=RunInfo.RegInfo.Cell.area;
centroid=RunInfo.RegInfo.Cell.centroid;

si_b=size(tcnorm_blank);
blank_resp = reshape(nanmean(tcnorm_blank((RunInfo.Noff+1):(RunInfo.Noff+4),:,:)),si_b(2),si_b(3));
%blank_resp = reshape((tcnorm_blank((RunInfo.Noff+RunInfo.Non+RunInfo.fbaseWindow),:,:)),si_b(2).*length(RunInfo.fbaseWindow),si_b(3));
blank_resp_SD = std(blank_resp);

dFF_sn_SD=RunInfo.Trace.dFFnoise;
nPhoton=RunInfo.Trace.nPhoton;
pt_shot_error=RunInfo.Trace.pt_shot_error;


for i=1:Ncells
    
   %Cell nanmean base luminance values
    OriStat.BasicStat(i).F = squeeze(nanmean(nanmean(base_data(:,:,:,i),2),3))-RunInfo.Trace.np_cont_ratio.*squeeze(nanmean(nanmean(base_data_np(:,:,:,i),2),3)); 
    
    OriStat.BasicStat(i).cell_contrast=squeeze(nanmean(nanmean(base_data(:,:,:,i),2),3))./squeeze(nanmean(nanmean(base_data_np(:,:,:,i),2),3));
   %response significance
   resp = squeeze(resp_data(:,:,:,i))-1;
   base = squeeze(base_data(:,:,:,i))-1;
   
   [OriStat.BasicStat(i).resp_p,h,OriStat.BasicStat(i).resp_stats]...
       = ranksum(base(:),resp(:));
   
   OriStat.BasicStat(i).RR=max(squeeze(nanmean(norm_data(:,:,:,i)-1,2)))./blank_resp_SD(i);
   
   OriStat.BasicStat(i).blank_resp_SD=blank_resp_SD(i);
   
   OriStat.BasicStat(i).pt_blank_resp_SD=blank_resp_SD(i)./sqrt(nreps);
   
   OriStat.BasicStat(i).dFF_sn_SD=dFF_sn_SD(i);
   OriStat.BasicStat(i).nPhoton=nPhoton(i);
   OriStat.BasicStat(i).pt_shot_error=pt_shot_error(i);
   
   OriStat.BasicStat(i).area=area(i);
   OriStat.BasicStat(i).x_cent=centroid(1,i);
   OriStat.BasicStat(i).y_cent=centroid(2,i);
   if size(centroid,1)==3
    OriStat.BasicStat(i).z_cent=centroid(3,i);
   end

   dir_table=OriStat.Ori(i).norm_table-1;   
   si=size(dir_table);
   ori_table=reshape(dir_table,si(1),si(2)./2,2);
   ori_table=squeeze(nanmean(ori_table,3));   
   
   %calculate ori_vector_tune_CI
   [LB68 LB95 UB68 UB95 CI68 CI95 boot_vector_tune]=ori_vector_tune_CI(ori_table);
   
   OriStat.Ori(i).ori_tune_bounds68=[LB68 UB68];
   OriStat.Ori(i).ori_tune_bounds95=[LB95 UB95];
   OriStat.Ori(i).ori_tune_CI68=CI68;
   OriStat.Ori(i).ori_tune_CI95=CI95;
   OriStat.Ori(i).ori_tune_boot=boot_vector_tune;
   
   OriStat.Ori(i).data_table_np=Ori_np(1).data_table;
   OriStat.Ori(i).norm_table_np=Ori_np(1).norm_table;
   OriStat.Ori(i).dir_ratio_change_np=Ori_np(1).dir_ratio_change;
   OriStat.Ori(i).ori_ratio_change_np=Ori_np(1).ori_ratio_change;
   
   OriStat.Ori(i).dir_vector_angle_np=Ori_np(1).dir_vector_angle;
   OriStat.Ori(i).dir_vector_mag_np=Ori_np(1).dir_vector_mag;
   OriStat.Ori(i).dir_vector_tune_np=Ori_np(1).dir_vector_tune;
   OriStat.Ori(i).dir_vector_tune_norm=OriStat.Ori(i).dir_vector_tune./Ori_np(1).dir_vector_tune;

   OriStat.Ori(i).ori_vector_angle_np=Ori_np(1).ori_vector_angle;
   OriStat.Ori(i).ori_vector_mag_np=Ori_np(1).ori_vector_mag;
   OriStat.Ori(i).ori_vector_tune_np=Ori_np(1).ori_vector_tune;  
   OriStat.Ori(i).ori_vector_tune_norm=OriStat.Ori(i).ori_vector_tune./Ori_np(1).ori_vector_tune;
   
    d= mod((OriStat.Ori(i).dir_vector_angle - OriStat.Ori(i).dir_vector_angle_np),360);
    if d > 180
        d=360-d;
    end
    OriStat.Ori(i).dir_vector_angle_npdiff=d;
    
    d= mod((OriStat.Ori(i).ori_vector_angle - OriStat.Ori(i).ori_vector_angle_np),180);
    if d > 90
        d=180-d;
    end
    OriStat.Ori(i).ori_vector_angle_npdiff=d; 
    
    OriStat.Ori(i).tcplot=squeeze(tcplot_tcnorm_ori(:,i));
end

OriStat.RunInfo=RunInfo;

OriStat.Code.doOriStatAK=getCode('doOriStatAK');

alloutfname=RunInfo.RegInfo.All.alloutfname;
save([alloutfname,'OriStat.mat'],'OriStat');
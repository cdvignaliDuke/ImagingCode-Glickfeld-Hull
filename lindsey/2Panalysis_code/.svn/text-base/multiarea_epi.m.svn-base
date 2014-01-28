%% df/f vs baseline f
base = 'G:\users\lindsey\analysisLG\active mice\';
date = '110515';
mouse = 'Y2';
data = '_run8_onlyon';
FOV = '_run9_FOV';

fn_data = fullfile(base, mouse, date,[date '_' mouse data '.tif']);
fn_FOV = fullfile(base, mouse, date,[date '_' mouse FOV '.tif']);

data_stack = readtiff(fn_data);
FOV_image = readtiff(fn_FOV);

LM_mean_data =(mean(mean(data_stack(221:240,271:290,:),2),1)-1)*100;
LM_FOV = mean(mean(FOV_image(221:240,271:290,:),2),1);
AL_mean_data = (mean(mean(data_stack(204:223,185:204,:),2),1)-1)*100;
AL_FOV = mean(mean(FOV_image(204:223,185:204,:),2),1);
RL_mean_data = (mean(mean(data_stack(179:198,114:133,:),2),1)-1)*100;
RL_FOV = mean(mean(FOV_image(179:198,114:133,:),2),1);
AM_mean_data = (mean(mean(data_stack(56:75,117:136,:),2),1)-1)*100;
AM_FOV = mean(mean(FOV_image(56:75,117:136,:),2),1);
PM_mean_data = (mean(mean(data_stack(34:53,153:172,:),2),1)-1)*100;
PM_FOV = mean(mean(FOV_image(34:53,153:172,:),2),1);

FOV = [LM_FOV AL_FOV RL_FOV AM_FOV PM_FOV]/max([LM_FOV AL_FOV RL_FOV AM_FOV PM_FOV]);
mean_data = [squeeze(LM_mean_data) squeeze(AL_mean_data) squeeze(RL_mean_data) squeeze(AM_mean_data) squeeze(PM_mean_data)];

figure;
for icond = 1:4
    subplot(2,2,icond);
    scatter(FOV, mean_data(icond,:))
    hold on;
    p = polyfit(FOV,mean_data(icond,:),1);
    x = [.5:.1:1];
    y = p(1)*x+p(2);
    plot(x,y);    
    ylim([0 10]);
end

%%  df/f vs position
base = 'G:\users\lindsey\analysisLG\active mice\';
date = '110518';
mouse = 'Y1';
data = '_5pos';

fn_data = fullfile(base, mouse, date,[date '_' mouse data '.tif']);
data_stack = readtiff(fn_data);

LM_mean_data =(mean(mean(data_stack(222:241, 301:320,:),2),1)-1)*100;
AL_mean_data =(mean(mean(data_stack(230:249,245:264,:),2),1)-1)*100;
PM_mean_data =(mean(mean(data_stack(16:35,189:208,:),2),1)-1)*100;

pos = [-500 -200 0 200 500];
mean_data = [squeeze(LM_mean_data)/max(squeeze(LM_mean_data)) squeeze(AL_mean_data)/max(squeeze(AL_mean_data)) squeeze(PM_mean_data)/max(squeeze(PM_mean_data))];
figure;
plot(pos,mean_data);
ylim([0 1])
box off


clear all
close all


load('/home/jason/data/iontoncm/all/B331_Z1260_ss001_rateanal.mat')

pre_rate1 = normrate;

pre_spon1 = sponrates;
pre_allspon1=allsponrates;

load('/home/jason/data/iontoncm/all/B331_Z1261_ss001_rateanal.mat')

ion_rate1 = normrate;
ion_spon1 =sponrates;
ion_allspon1=allsponrates;

hold on
plot(mean(pre_rate1))
plot(mean(ion_rate1),'r')

plot(mean(pre_rate1)+(std(pre_rate1)./(size(pre_rate1,1))),'--' )
plot(mean(pre_rate1)-(std(pre_rate1)./(size(pre_rate1,1))),'--' )

plot(mean(ion_rate1)+(std(ion_rate1)./(size(ion_rate1,1))),'--r' )
plot(mean(ion_rate1)-(std(ion_rate1)./(size(ion_rate1,1))),'--r' )
xlabel('Song')
ylabel('Firing rate')
legend({'Control'; 'Gabazine'})
title('Cell 1');
hold off



load('/home/jason/data/iontoncm/all/B331_Z3020_ss001_rateanal.mat')

pre_rate2 = normrate;


load('/home/jason/data/iontoncm/all/B331_Z3021_ss001_rateanal.mat')


ion_rate2 = normrate;
figure
hold on
plot(mean(pre_rate2))
plot(mean(ion_rate2),'r')

plot(mean(pre_rate2)+(std(pre_rate2)./(size(pre_rate2,1))),'--' )
plot(mean(pre_rate2)-(std(pre_rate2)./(size(pre_rate2,1))),'--' )

plot(mean(ion_rate2)+(std(ion_rate2)./(size(ion_rate2,1))),'--r' )
plot(mean(ion_rate2)-(std(ion_rate2)./(size(ion_rate2,1))),'--r' )
xlabel('Song')
ylabel('Firing rate')
legend({'Control'; 'Gabazine'})
title('Cell 2');

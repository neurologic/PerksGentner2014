function plotiontoallsongsumm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plots mean firing rates to each song stimulus for baseline, gabazine and
%recovery conditions. Makes seperate plots for normalized and standard rate
file=dir('*iontodata.mat');
load(file.name);


figure
plot(baserate,'-ok')
hold on
plot(gzrate,'-or')
plot(recrate,'-og')
ylabel('Firing Rate (spikes/sec)')
set(gca,'xticklabel',stims)
xlabel('Song')

figure
[x,y]=max(baserate);
normbaserate=100*baserate./(baserate(y));
normgzrate=100*gzrate./(gzrate(y));
normrecrate=100*recrate./(recrate(y));
plot(normbaserate,'-ok')
hold on
plot(normgzrate,'-or')
plot(normrecrate,'-og')
ylabel('% of Maximum Firing Rate')
set(gca,'xticklabel',stims)
xlabel('Song')


end
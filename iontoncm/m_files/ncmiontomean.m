cd ~/data/iontoncm/mean
clear all 
close all

files=dir('*ncmionto.mat');

for i =1:length(files)
   load(files(i).name);
   
   mbasebias(i)=basebias;
   mdrugbias(i)=drugbias;
   mwashbias(i)=washbias;
   
   mbasesparse(i)=baselifesparse;
   mdrugsparse(i)=druglifesparse;
   mwashsparse(i)=washlifesparse;

   meanlearnperinc(i)=mean(l_songperinc);
   meanunfamperinc(i)=mean(u_songperinc);
   
   mbasego(i)=zbasego;
   mbasenogo(i)=zbasenogo;
   mbasenov(i)=zbasenov;
   mbasetrain(i)=zbasetrain;
   
   mdruggo(i)=zdruggo;
   mdrugnogo(i)=zdrugnogo;
   mdrugnov(i)=zdrugnov;
   mdrugtrain(i)=zdrugtrain;
   
    depths(i)=depth;
    
    mtraindiff(i)=mean(trainperdiff);
    mnovdiff(i)=mean(novperdiff);
    
    mbaseisi(i)=mean(baseisis*1000); % in msec
    mdrugisi(i)=mean(drugisis*1000);
    mwashisi(i)=mean(washisis*1000);
    
end

figure
   bar([1 2 3],[mean(mbasego) mean(mbasenogo) mean(mbasenov)],'b')
hold on
   errorbar([1 2 3],[mean(mbasego) mean(mbasenogo) mean(mbasenov)],[sem(mbasego) sem(mbasenogo) sem(mbasenov)],'.b')
  bar([5 6 7],[mean(mdruggo) mean(mdrugnogo) mean(mdrugnov)],'r')
      errorbar([5 6 7],[mean(mdruggo) mean(mdrugnogo) mean(mdrugnov)],[sem(mdruggo) sem(mdrugnogo) sem(mdrugnov)],'.r')
   xlim([0 8])
   set(gca,'xtick',[1 2 3 4 5 6 7 ])
  set(gca,'xticklabel',{'base Go'; 'base Nogo'; 'base Nov'; ' '; 'drug Go'; 'drug Nogo'; 'drug Nov'})
  ylabel('zscore firing rate')
hold off

figure
   bar([1 2],[mean(mbasetrain) mean(mbasenov)],'b')
hold on
   errorbar([1 2],[mean(mbasetrain) mean(mbasenov)],[sem(mbasetrain) sem(mbasenov)],'.b')
  bar([4 5] ,[mean(mdrugtrain) mean(mdrugnov)],'r')
      errorbar([4 5],[mean(mdrugtrain) mean(mdrugnov)],[sem(mdrugtrain) sem(mdrugnov)],'.r')
   xlim([0 6])
   set(gca,'xtick',[1 2 3 4 5 ])
  set(gca,'xticklabel',{'base train'; 'base Nov'; ' '; 'drug train'; 'drug Nov'})
  ylabel('zscore firing rate')
hold off
   








figure
plot([mbasesparse; mdrugsparse; mwashsparse],'-o')
ylabel('lifetime sparseness')
xlim([0.5 3.5 ])
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',{'baseline';'drug'; 'wash'})

figure
   bar([mean(mtraindiff) mean(mnovdiff)],'k')
hold on
   errorbar([mean(mtraindiff) mean(mnovdiff)],[sem(mtraindiff) sem(mnovdiff)],'.k')
  xlim([0 3])
  set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
  ylabel('% change in firing rate with gabazine')
  
  hold off

  
  figure
   bar([mean(mbaseisi) mean(mdrugisi)],'k')
   hold on
   errorbar([mean(mbaseisi) mean(mdrugisi)],[sem(mbaseisi) sem(mdrugisi)],'.k')

ylabel('ISI (msec)')
xlim([0 3 ])
%set(gca,'xtick',[1 2 ])
set(gca,'xticklabel',{'baseline';'drug'})


figure
plot([mbasebias; mdrugbias; mwashbias],'-o')
ylabel('unfam-learned bias')
xlim([0.5 3.5])
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',{'baseline';'drug'; 'wash'})

jo=find(mbasebias<0);

nmbasebias=mbasebias;
nmbasebias(jo)=abs(nmbasebias(jo));

nmdrugbias=mdrugbias;
nmdrugbias(jo)=abs(nmdrugbias(jo));

figure
plot([abs(mbasebias); abs(mdrugbias); abs(mwashbias)],'-o')
ylabel('absolute unfam-learned bias')
xlim([0.5 3.5])
ylim([0 .6])
set(gca,'xtick',[1 2 3])
set(gca,'xticklabel',{'baseline';'drug'; 'wash'})

% %only use sites that prefer novel songs in baseline condition
% po=find(mbasebias>0);
% 
% figure
% plot([mbasebias(po); mdrugbias(po); mwashbias(po)],'-o')
% ylabel('unfam-learned bias')
% xlim([0.5 3.5])
% set(gca,'xtick',[1 2 3])
% set(gca,'xticklabel',{'baseline';'drug'; 'wash'})
% title('sites that prefer unfamiliar songs')
% 
% 
% %sites that prefer learned songs
% mo=find(mbasebias<0);
% figure
% plot([mbasebias(mo); mdrugbias(mo); mwashbias(mo)],'-o')
% ylabel('unfam-learned bias')
% xlim([0.5 3.5])
% set(gca,'xtick',[1 2 3])
% set(gca,'xticklabel',{'baseline';'drug'; 'wash'})
% title('sites that prefer learned songs')
% 











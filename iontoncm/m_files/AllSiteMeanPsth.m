% Makes Average PSTH for all sites where 5 sec stims were used

cd ~/data/iontoncm/splitbydrug/all/
files=dir('*iontodata.mat');

Ivals=[];
SiteInd=1;
for i=1:length(files)
   load(files(i).name)
   stimlengths{i}=stimlength;
   if strmatch(stimlength,'s')==1
      AllSongPsth.Base(SiteInd,:) = mean(Psth.Base);
      AllSongPsth.Gz(SiteInd,:) = mean(Psth.Gz);
      AllSongPsth.Rec(SiteInd,:) = mean(Psth.Rec);
      
      LearnPsth.Base(SiteInd,:) = mean(Psth.Base(trainstims,:));
      LearnPsth.Gz(SiteInd,:) = mean(Psth.Gz(trainstims,:));
      LearnPsth.Rec(SiteInd,:) = mean(Psth.Rec(trainstims,:));
      
      UnfamPsth.Base(SiteInd,:) = mean(Psth.Base(novstims,:));
      UnfamPsth.Gz(SiteInd,:) = mean(Psth.Gz(novstims,:));
      UnfamPsth.Rec(SiteInd,:) = mean(Psth.Rec(novstims,:));
      
      MaxFR = max(max(Psth.Base));
      NormAllSongPsth.Base(SiteInd,:) = mean(Psth.Base)./MaxFR;
      NormAllSongPsth.Gz(SiteInd,:) = mean(Psth.Gz)./MaxFR;
      NormAllSongPsth.Rec(SiteInd,:) = mean(Psth.Rec)./MaxFR;
      
      NormLearnPsth.Base(SiteInd,:) = mean(Psth.Base(trainstims,:))./MaxFR;
      NormLearnPsth.Gz(SiteInd,:) = mean(Psth.Gz(trainstims,:))./MaxFR;
      NormLearnPsth.Rec(SiteInd,:) = mean(Psth.Rec(trainstims,:))./MaxFR;
      
      NormUnfamPsth.Base(SiteInd,:) = mean(Psth.Base(novstims,:))./MaxFR;
      NormUnfamPsth.Gz(SiteInd,:) = mean(Psth.Gz(novstims,:))./MaxFR;
      NormUnfamPsth.Rec(SiteInd,:) = mean(Psth.Rec(novstims,:))./MaxFR;
      
      Ivals=[Ivals i];
      Depths(SiteInd) = depth;
      for j = 1:size(Psth.Base,1)
        [r,p]=corrcoef(Psth.Base(j,:),Psth.Gz(j,:));
        Rvals(j)=r(2);
      end
      RGzBase(SiteInd)=nanmean(Rvals);
      
       SiteInd=SiteInd+1;
       
   end





end

Dorsal = find(Depths < 2580);
Ventral = find(Depths >= 2580);
%%
figure
subplot(2,3,1);
plot(xbins,mean(AllSongPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(AllSongPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec All Songs')
title('Dorsal')
hold off
subplot(2,3,4);
plot(xbins,mean(AllSongPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(AllSongPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec All Songs')
title('Ventral')
hold off
subplot(2,3,2);
plot(xbins,mean(LearnPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(LearnPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec Learned Songs')
title('Dorsal')
hold off
subplot(2,3,3);
plot(xbins,mean(UnfamPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(UnfamPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec Unfamiliar Songs')
title('Dorsal')
hold off
subplot(2,3,5);
plot(xbins,mean(LearnPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(LearnPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec Learned Songs')
title('Ventral')
hold off
subplot(2,3,6);
plot(xbins,mean(UnfamPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(UnfamPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Mean Spikes/Sec Unfamiliar Songs')
title('Ventral')
hold off

%%
figure
subplot(2,3,1);
plot(xbins,mean(NormAllSongPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(NormAllSongPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec All Songs')
title('Dorsal')
hold off
subplot(2,3,4);
plot(xbins,mean(NormAllSongPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(NormAllSongPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec All Songs')
title('Ventral')
hold off
subplot(2,3,2);
plot(xbins,mean(NormLearnPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(NormLearnPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec Learned Songs')
title('Dorsal')
hold off
subplot(2,3,3);
plot(xbins,mean(NormUnfamPsth.Base(Dorsal,:)),'k')
hold on
plot(xbins,mean(NormUnfamPsth.Gz(Dorsal,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec Unfamiliar Songs')
title('Dorsal')
hold off
subplot(2,3,5);
plot(xbins,mean(NormLearnPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(NormLearnPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec Learned Songs')
title('Ventral')
hold off
subplot(2,3,6);
plot(xbins,mean(NormUnfamPsth.Base(Ventral,:)),'k')
hold on
plot(xbins,mean(NormUnfamPsth.Gz(Ventral,:)),'r')
xlabel('Time');
ylabel('Norm Spikes/Sec Unfamiliar Songs')
title('Ventral')
hold off

%%
figure 
subplot(1,2,1);
plot(xbins,mean(NormLearnPsth.Gz(Dorsal,:)-NormLearnPsth.Base(Dorsal,:)),'m') 
hold on
plot(xbins,mean(NormUnfamPsth.Gz(Dorsal,:)-NormUnfamPsth.Base(Dorsal,:)),'b') 
xlabel('Time');
ylabel('Difference Gabazine - Baseline Spikes/Sec')
title('Dorsal')
legend({'Learned';'Unfamiliar'});
hold off
subplot(1,2,2);
plot(xbins,mean(NormLearnPsth.Gz(Ventral,:)-NormLearnPsth.Base(Ventral,:)),'m') 
hold on
plot(xbins,mean(NormUnfamPsth.Gz(Ventral,:)-NormUnfamPsth.Base(Ventral,:)),'b') 
xlabel('Time');
ylabel('Difference Gabazine - Baseline Spikes/Sec')
title('Ventral')
legend({'Learned';'Unfamiliar'});
hold off

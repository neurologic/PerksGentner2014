function MakeMotifFRColorMapIonto(infile)

load(infile)
numFRbins = 7;
clear sort;

minFRA = min(min([basemotreprate gzmotreprate recmotreprate]));
maxFRA = max(max([basemotreprate gzmotreprate recmotreprate]));
minFRBase = min(min([basemotreprate ]));
maxFRBase = max(max([basemotreprate ]));
minFRGz = min(min([ gzmotreprate ]));
maxFRGz = max(max([ gzmotreprate ]));
minFRRec = min(min([ recmotreprate]));
maxFRRec = max(max([ recmotreprate]));




binsA = linspace(minFRA, maxFRA, numFRbins);
binsBase=linspace(minFRBase,maxFRBase,numFRbins);
binsGz=linspace(minFRGz,maxFRGz,numFRbins);
binsRec=linspace(minFRRec,maxFRRec,numFRbins);

baseMotHist=zeros(size(basemotreprate,2),length(binsA));
baseMotHistBaseBins=zeros(size(baseMotHist));
basemotreprate=basemotreprate';
for i=1:size(basemotreprate,1)
   baseMotHist(i,:)=histc(basemotreprate(i,:),binsA);
   baseMotHistBaseBins(i,:)=histc(basemotreprate(i,:),binsBase);
           
end

gzMotHist=zeros(size(gzmotreprate,2),length(binsA));
gzmotreprate=gzmotreprate';
gzMotHistGzBins=zeros(size(gzMotHist));
for i=1:size(gzmotreprate,1)
   gzMotHist(i,:)=histc(gzmotreprate(i,:),binsA);
   gzMotHistGzBins(i,:)=histc(gzmotreprate(i,:),binsGz);    
end

recMotHist=zeros(size(recmotreprate,2),length(binsA));
recMotHistRecBins=zeros(size(recMotHist));
recmotreprate=recmotreprate';
for i=1:size(recmotreprate,1)
   recMotHist(i,:)=histc(recmotreprate(i,:),binsA);
   recMotHistRecBins(i,:)=histc(recmotreprate(i,:),binsRec);   
end
sumMat=[];
for j=1:size(baseMotHist,2)
    
    
    sumMat=[sumMat sum(baseMotHist,2)] ;

end

baseMotHist=baseMotHist./sumMat;
gzMotHist=gzMotHist./sumMat;
recMotHist=recMotHist./sumMat;

baseMotHistBaseBins=baseMotHistBaseBins./sumMat;
gzMotHistGzBins=gzMotHistGzBins./sumMat;
recMotHistRecBins=recMotHistRecBins./sumMat;


[C,I]=sort(baseMotHist(:,1));
baseMotHist=baseMotHist(flipud(I),:);
[C,I]=sort(gzMotHist(:,1));
gzMotHist=gzMotHist(flipud(I),:);
[C,I]=sort(recMotHist(:,1));
recMotHist=recMotHist(flipud(I),:);

[C,I]=sort(baseMotHistBaseBins(:,1));
baseMotHistBaseBins=baseMotHistBaseBins(flipud(I),:);
[C,I]=sort(gzMotHistGzBins(:,1));
gzMotHistGzBins=gzMotHistGzBins(flipud(I),:);
[C,I]=sort(recMotHistRecBins(:,1));
recMotHistRecBins=recMotHistRecBins(flipud(I),:);

figure 
subplot(1,3,1); imagesc(baseMotHist)
colormap(hot); xlabel('Spikes/Sec'); ylabel('Motif Number');title('Baseline');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsA));
subplot(1,3,2); imagesc(gzMotHist)
colormap(hot);xlabel('Spikes/Sec'); ylabel('Motif Number');title('Gabazine');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsA));
subplot(1,3,3); imagesc(recMotHist)
colormap(hot);xlabel('Spikes/Sec'); ylabel('Motif Number');title('Recovery');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsA));
%colorbar

figure 
subplot(1,3,1); imagesc(baseMotHistBaseBins)
colormap(hot); xlabel('Spikes/Sec'); ylabel('Motif Number');title('Baseline');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsBase));
subplot(1,3,2); imagesc(gzMotHistGzBins)
colormap(hot);xlabel('Spikes/Sec'); ylabel('Motif Number');title('Gabazine');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsGz));
subplot(1,3,3); imagesc(recMotHistRecBins)
colormap(hot);xlabel('Spikes/Sec'); ylabel('Motif Number');title('Recovery');
set(gca,'xtick',[ 1 2 3 4 5 6 7]); set(gca,'xticklabel', round(binsRec));
%colorbar

figure
[X,Y]=sort(basemotrate,'descend');
subplot(1,3,1); plot(basemotrate(Y),'-ok')
hold on
plot(gzmotrate(Y),'-or')
plot(recmotrate(Y),'-og')
ylabel('Firing Rate (Spikes/Sec)')
xlabel('Motif')
legend({'Baseline';'Gabazine';'Recovery'})

[X,Y]=sort(basemotrate,'descend');
subplot(1,3,2); plot(basemotrate(Y),'-ok')
hold on
[X,Y]=sort(gzmotrate,'descend');
plot(gzmotrate(Y),'-or')
[X,Y]=sort(recmotrate,'descend');
plot(recmotrate(Y),'-og')
ylabel('Firing Rate (Spikes/Sec)')
xlabel('Motif')
legend({'Baseline';'Gabazine';'Recovery'})

subplot(1,3,3); [X,Y]=sort(basemotrate,'descend');
plot(100*(basemotrate(Y)./max(basemotrate)),'-ok')
hold on
[X,Y]=sort(gzmotrate,'descend');
plot(100*(gzmotrate(Y)./max(gzmotrate)),'-or')
[X,Y]=sort(recmotrate,'descend');
plot(100*(recmotrate(Y)./max(recmotrate)),'-og')
ylabel('% of Maximum Firing Rate ')
xlabel('Motif')
legend({'Baseline';'Gabazine';'Recovery'})
end













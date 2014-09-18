load PMDdata % PMD dataset from monkey G
%%
icell = 1

figure;

for itrial = 1:size(PMDdata(1).spikes,1)
    thesespikes = find(PMDdata(icell).spikes(itrial,:));
    for ispike = 1:size(thesespikes,2)
        
        line([thesespikes(ispike) thesespikes(ispike)], [double(itrial)-.4 double(itrial)+.4],'Linewidth',.5,'Color','k');
                                                               
    end     
   
end
%%
hfig = figure;
 cond = 200;
data = base_data_All;
for itrial = 1:size(data(cond).spikes,1)
    thesespikes = find(data(cond).spikes(itrial,:));
  
    for ispike = 1:size(thesespikes,2)
        
        line([thesespikes(ispike) thesespikes(ispike)], [double(itrial)-.4 double(itrial)+.4],'Linewidth',.5,'Color','k');
                                                                    
    end     
   
end

set(hfig,'Position',pos);
%%


files=dir('*dbmean_toe.txt');
fileind = 1;
for i=1:length(files)
    close all
    clearvars -except PMDdata data base_data gz_data dosave files i start stop blocksize blocktype sitetype smooth doprintout
    
%  batchrasterionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
 [base, gz]= FanoFormat_ionto(files(i).name,-2,12,10,[1 2 5],'mu',0,0);
 
 base_data(fileind).spikes = base;
 gz_data(fileind).spikes = gz;
 fileind = fileind + 1;
end

% need to split the data up so that I get baseline trials for one data struct and
% gabazine trials for the other data struct

%%
figure;
for itrial = 1:size(base_data(1),2)
    thesespikes = cell2mat(toes{itrial});
     for ispike = 1:size(thesespikes,1)  
        line([thesespikes(ispike) thesespikes(ispike)], [double(itrial)-.4 double(itrial)+.4],'Linewidth',.5,'Color','k');                                      
    end     
end
%%
 % target onset is at 400 ms
%alldelays>= 400ms
times = 50:50:6000; % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 2000; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
fanoParams.matchReps = 0;
fanoParams.binSpacing = 0.25;
fanoParams.binWidth = 50;
Result_base = VarVsMean(base_data, times, fanoParams); 
plotFano(Result_base);

Result_gz = VarVsMean(gz_data, times, fanoParams); 
plotFano(Result_gz);

figure;hold on
plot(Result_base.meanRateAll,'color','k')
plot(Result_gz.meanRateAll,'color','r')
axis tight


figure;hold on
plot(Result_base.FanoFactorAll,'color','k')
plot(Result_gz.FanoFactorAll,'color','r')
axis tight
    
%%
 % target onset is at 400 ms
%alldelays>= 400ms
times = 200:25:850; % from 200 ms before target onset until 450 ms after.
fanoParams.alignTime = 400; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
Result = VarVsMean(PMDdata, times, fanoParams); plotFano(Result);
%%
scatterParams.axLim = 'auto';
scatterParams.axLen = 5;
plotScatter(Result, -100, scatterParams);
text(2.5, 7, '100 ms before target', 'hori', 'center'); 
plotScatter(Result, 100, scatterParams);
text(2.5, 7, '100 ms after target', '?hori', '?center'); 
plotScatter(Result, 300, scatterParams);
text(2.5, 7, '300 ms after target', '?hori', '?center');
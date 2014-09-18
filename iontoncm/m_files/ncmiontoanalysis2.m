%%calculate the selectivity index for a single cell
% For cells with no silence response file
clear all
close all

!rm *_ncmionto.mat

doprintout=0;
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'};

bird = {'st444'; 'st445'; 'st493'};
trainset = {'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'};
goset = {'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'};
nogoset = { 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'};
novset = { 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'};


basereps=5;
washreps=5;
drugreps=5;


pretime=2; %in sec
posttime=2; % in sec
binsize = 50; %in msec

files = dir('*dbmean_toe.txt');

for i=1:length(files)
    toelist = files(i).name;
    %toelist = char(files(i));
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);  
    nreps=nreps;
end
  subindx=strcmpi(bird,subjectID);
 
bps = 1000/binsize;

%oldrate=zeros(nreps,length(files));

allspikes=[];
trainstims=[];
novstims=[];
gostims=[];
nogostims=[];
for i=1:length(files)
    toelist = files(i).name;
    %toelist = char(files(i));
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);  
    stim = strtok(stimfile, '_');   
    stim = strtok(stim, '.');
     if sum(strcmpi(stim,trainset(subindx,:)))>0
    trainstims=[trainstims  i];
    end
    
    if sum(strcmpi(stim,goset(subindx,:)))>0
    gostims=[gostims  i];
    end
    
    if sum(strcmpi(stim,nogoset(subindx,:)))>0
    nogostims=[nogostims  i];
    end
    if sum(strcmpi(stim,novset(subindx,:)))>0
    novstims=[novstims  i];
    end
        
    stimdur = stimtimes(strmatch(stim, stimnames, 'exact'));
    stimdurs(i)=stimdur;
    othertoes = alltoes(find(alltoes>=0));
    othertoes = othertoes(find(othertoes<=stimdur));
    xbins = linspace(-2, stimdur+2, (stimdur+4)*bps);
    allxbins{i}=xbins;
    numxbins(i)=length(xbins);
    %psth{i} = histc(othertoes, xbins);
    psth=[];
    
    for moo=1:nreps
            
            newtoes = toes{moo}{1}(find(toes{moo}{:}>=0));%only use spikes after stim starts
            newtoes = newtoes(find(newtoes<=stimdur)); %get rid of spikes after stim over
            %newtoes = round(newtoes*1000);  %may want to take this out
            %newtoes = (newtoes*1000); 
            %tpsth=zeros(size(xbins));
            atoes=toes{moo}{1};
            
            if isempty(atoes)==1
                
                tpsth=zeros(1,length(xbins));
            else 
                
                tpsth=histc(atoes, xbins);
            
                
                
            end
            psth(moo,:) = tpsth';
            
            oldrate(moo,i) = length(newtoes)/(stimdur);
                
               
              sponr(moo,i) = length( toes{moo}{1}(find(toes{moo}{:}<0)))/2;
              
         
    end
    
    
    allpsth{i}=psth;
    numspikes(i)=length(alltoes);
    spontoes = alltoes(find(alltoes<0));
    sponrate = length(spontoes)/(2*double(nreps)); 
   
   % normrate = oldrate/sponrate;
    normrate = oldrate;
    stims{i} = stim;
    reps(i) = nreps;
    %normrates(i) = normrate;
    sponrates(i)=sponrate;

end

purerate=mean(normrate);
muspon = mean(sponrates);
threshrate=mean(sponrates)+std(sponrates);
numsigstims=length(find(mean(normrate)>threshrate));
psigstims=numsigstims./length(stims);

tonefiles = dir('*khz_toe.txt');

for k=1:length(tonefiles)
    toelist = tonefiles(k).name;
    %toelist = char(files(i));
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);  
    stim = strtok(stimfile, '_');   
    stim = strtok(stim, '.');
        
    tonestimdur = stimtimes(strmatch(stim, stimnames, 'exact'));
    othertoes = alltoes(find(alltoes>=0));
    othertoes = othertoes(find(othertoes<=tonestimdur));
    xbins = linspace(-2, tonestimdur+2, (tonestimdur+4)*bps);
    toneallxbins{k}=xbins;
    numxbins(k)=length(xbins);
    %psth{i} = histc(othertoes, xbins);
    tonepsth=[];
    
    for joo=1:nreps
            
            tnewtoes = toes{joo}{1}(find(toes{joo}{:}>=0));%only use spikes after stim starts
            tnewtoes = tnewtoes(find(tnewtoes<=tonestimdur)); %get rid of spikes after stim over
            %newtoes = round(newtoes*1000);  %may want to take this out
            %newtoes = (newtoes*1000); 
            %tpsth=zeros(size(xbins));
            atoes=toes{joo}{1};
            
            if isempty(atoes)==1
                
                tpsth=zeros(1,length(xbins));
            else 
                
                tpsth=histc(atoes, xbins);
            
                
                
            end
            tonepsth(joo,:) = tpsth';
            
            toneoldrate(joo,k) = length(tnewtoes)/(tonestimdur);
                
               
              tonesponr(joo,k) = length( toes{joo}{1}(find(toes{joo}{:}<0)))/2;
              
         
    end
    
    
    toneallpsth{k}=tonepsth;
    %numspikes(k)=length(alltoes);
    spontoes = alltoes(find(alltoes<0));
    sponrate = length(spontoes)/(2*double(nreps)); 
   
   % normrate = oldrate/sponrate;
    %normrate = oldrate;
    tonestims{k} = stim;
    tonereps(k) = nreps;
    %normrates(k) = normrate;
    tonesponrates(k)=sponrate;

end


l_songperinc=[];
u_songperinc=[];
allperinc=[];
for j=1:size(allpsth,2)
   currpsth=allpsth{j};
   
   basepsth=mean(currpsth(1:basereps,:));
   drugpsth=mean(currpsth(basereps+1:nreps-washreps,:));
   washpsth=mean(currpsth(nreps-washreps+1:nreps,:));
   mcurrpsth=mean(currpsth);
   
   %minval=min(mcurrpsth(find(mcurrpsth>0)));
   minval=0.1;
   
   figure
   hold on;
   plot(allxbins{j},basepsth)
   plot(allxbins{j},drugpsth,'r')   
   plot(allxbins{j},washpsth,'k') 
   xlabel('time (sec)');
   ylabel('mean spikes/sec');
   legend({'baseline';'drug';'washout'})
   figname = strrep(sprintf('subj:%s  pen:%s  site:%s  sort:%s  stim:%s', subjectID,pen,site,sort,stims{j}), '_', '\_');
   title(figname);
   hold off
   
   set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
   if doprintout==1
       
        print -dpsc2 -cmyk
       % print -noui 
    end
   
   
 
    drugpsth(find(drugpsth==0))=minval;
    basepsth(find(basepsth==0))=minval;
    perinc= 100*((drugpsth-basepsth)./(basepsth));
%   
   allperinc=[allperinc perinc]; 
   sponbins=find(allxbins{j}<=0);
   songbins=find(allxbins{j}>0&allxbins{j}<=stimdurs(j));
   
   sponperinc(j)=mean(perinc(sponbins));
   songperinc(j)=mean(perinc(songbins));
   
 if sum(find(j==trainstims))>0
   l_songperinc=[l_songperinc mean(perinc(songbins))];
 else
   u_songperinc=[u_songperinc mean(perinc(songbins))];
 end

   spondrugpsth=drugpsth(sponbins);
   songdrugpsth=drugpsth(songbins);
   
   sponbasepsth=basepsth(sponbins);
   songbasepsth=basepsth(songbins);
   
   sponperincarea(j)=100*(trapz(spondrugpsth)-trapz(sponbasepsth))./(trapz(sponbasepsth));
   songperincarea(j)=100*(trapz(songdrugpsth)-trapz(songbasepsth))./(trapz(songbasepsth));   
   
end


figure
hist(perinc,50)
xlabel('% Change from baseline firing rate')
ylabel('# of 50ms bins')
figname = strrep(sprintf('subj:%s  pen:%s  site:%s  sort:%s', subjectID,pen,site,sort), '_', '\_');
title(figname);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
if doprintout==1
    print -dpsc2 -cmyk
    %print -noui 
end

% figure
% hold on
% bar([mean(sponperinc) mean(songperinc)],'k')
% errorbar([mean(sponperinc) mean(songperinc)],[sem(sponperinc) sem(songperinc)],'.k')
% ylabel('% change from baseline firing rate')
% set(gca,'xtick',[1 2])
% set(gca,'xticklabel',{'Spontaneous'; 'Driven'});
% xlim([0 3])
% figname = strrep(sprintf('subj:%s  pen:%s  site:%s  sort:%s', subjectID,pen,site,sort), '_', '\_');
% title(figname);
% hold off
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperOrientation', 'landscape');
% set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
% if doprintout==1
%     print -noui 
% end

figure
hold on
bar([mean(sponperincarea) mean(songperincarea)],'k')
errorbar([mean(sponperincarea) mean(songperincarea)],[sem(sponperincarea) sem(songperincarea)],'.k')
ylabel('% change in area under psth compared to baseline')
set(gca,'xtick',[1 2])
set(gca,'xticklabel',{'Spontaneous'; 'Driven'});
xlim([0 3])
figname = strrep(sprintf('subj:%s  pen:%s  site:%s  sort:%s', subjectID,pen,site,sort), '_', '\_');
title(figname);
hold off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
if doprintout==1
    print -dpsc2 -cmyk
    %print -noui 
end


if isempty(strfind(pen,'Right'))==0||isempty(strfind(pen,'RIGHT'))==0
    
   hemi='R';
elseif isempty(strfind(pen,'Left'))==0||isempty(strfind(pen,'LEFT'))==0
    
    hemi='L';
end
 
ymax=round(max(oldrate(:)));





figure
% bar([11 40],[10 80], [.5 .5 .5])
hold on
patch([basereps basereps nreps-washreps  nreps-washreps],[0 ymax ymax 0],[.5 .5 .5],'edgealpha',0)
plot(oldrate)
ylim([0 ymax])
xlabel('repetition')
ylabel('mean firing rate (spikes/sec)')
figname = strrep(sprintf('subj:%s  pen:%s  site:%s  sort:%s', subjectID,pen,site,sort), '_', '\_');
title(figname);
hold off
set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
   if doprintout==1
       %print -dpsc2 -rgb
        print -noui 
   end

allbasenov=mean(oldrate(1:basereps,novstims));
basenov=mean(allbasenov);
alldrugnov=mean(oldrate(basereps+1:basereps+drugreps,novstims));
drugnov=mean(alldrugnov);
allwashnov=mean(oldrate(basereps+drugreps+1:basereps+drugreps+washreps,novstims));
washnov=mean(allwashnov);
allbasetrain=mean(oldrate(1:basereps,trainstims));
basetrain=mean(allbasetrain);
alldrugtrain=mean(oldrate(basereps+1:basereps+drugreps,trainstims));
drugtrain=mean(alldrugtrain);
allwashtrain=mean(oldrate(basereps+drugreps+1:basereps+drugreps+washreps,trainstims));
washtrain=mean(allwashtrain);


basebias=(basenov-basetrain)./(basenov+basetrain)
drugbias=(drugnov-drugtrain)./(drugnov+drugtrain)
washbias=(washnov-washtrain)./(washnov+washtrain)


baserate=mean(oldrate(1:basereps,:));
drugrate=mean(oldrate(basereps+1:basereps+drugreps,:));
washrate=mean(oldrate(basereps+drugreps+1:basereps+drugreps+washreps,:));

% figure
% plot(baserate,'-o')
% hold on
% plot(drugrate,'-or')
% plot(washrate,'-ok')
% set(gca,'xticklabel',stims)
% legend({'baseline';'drug';'washout'})
% ylabel('firing rate (spikes/sec)')
% hold off
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperOrientation', 'landscape');
% set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
% if doprintout==1
%     print -noui 
% end  

figure
[x1,y1]=max(baserate);
nbaserate=100*(baserate./baserate(y1));
ndrugrate=100*(drugrate./drugrate(y1));
nwashrate=100*(washrate./washrate(y1));
plot(nbaserate,'-o')
hold on
plot(ndrugrate,'-or')
plot(nwashrate,'-ok')
set(gca,'xticklabel',stims)
legend({'baseline';'drug';'washout'})
ylabel('normalized firing rate');
hold off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
if doprintout==1
    %print -noui 
    print -dpsc2 -cmyk
end  



if isempty(tonefiles)==0
    tonebaserate=mean(toneoldrate(1:basereps,:));
    tonedrugrate=mean(toneoldrate(basereps+1:basereps+drugreps,:));
    tonewashrate=mean(toneoldrate(basereps+drugreps+1:basereps+drugreps+washreps,:));
    
    toneindex=[3 4 5 6 1 2]; %put tones in order
%     figure
%     plot(tonebaserate(toneindex),'-o')
%     hold on
%     plot(tonedrugrate(toneindex),'-or')
%     plot(tonewashrate(toneindex),'-ok')
%     set(gca,'xtick',1:6);
%     set(gca,'xticklabel',tonestims(toneindex))
%     legend({'baseline';'drug';'washout'})
%     hold off
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf, 'PaperUnits', 'inches');
%     set(gcf, 'PaperOrientation', 'landscape');
%     set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
%     xlabel('tone (khz)')
%     ylabel('firing rate (spikes/sec)')
%     if doprintout==1
%         print -noui 
%     end  
    
    
    figure
    [x2,y2]=max(tonebaserate);
    tonebaserate=(tonebaserate./tonebaserate(y2))*100;
    tonedrugrate=(tonedrugrate./tonedrugrate(y2))*100;
    tonewashrate=(tonewashrate./tonewashrate(y2))*100;
    plot(tonebaserate(toneindex),'-o')
    hold on
    plot(tonedrugrate(toneindex),'-or')
    plot(tonewashrate(toneindex),'-ok')
    set(gca,'xtick',1:6);
    set(gca,'xticklabel',tonestims(toneindex))
    legend({'baseline';'drug';'washout'})
    hold off
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperPosition', [.15 .15 10.5 8]);
    xlabel('tone (khz)')
    ylabel('normalized firing rate ')
    if doprintout==1
       % print -noui 
        print -dpsc2 -cmyk
    end  
end



basenum=(sum(baserate./(length(baserate)))).^2;
baseden=sum((baserate.^2)./(length(baserate)));
baselifesparse=1-(basenum/baseden);

drugnum=(sum(drugrate./(length(drugrate)))).^2;
drugden=sum((drugrate.^2)./(length(drugrate)));
druglifesparse=1-(drugnum/drugden);

washnum=(sum(washrate./(length(washrate)))).^2;
washden=sum((washrate.^2)./(length(washrate)));
washlifesparse=1-(washnum/washden)



    [xd,yd]=strtok(pen,'_');
[xd,zd]=strtok(yd,'_');
[RC,qd]=strtok(xd,'R');

[ML,nd]=strtok(zd,'_L');
    
    
    
    depth=str2double(strtok(site,'Z'));
% [T,R]=strtok(pwd, '/');
% [T,R]=strtok(R, '/');
% [T,R]=strtok(R, '/');
% [T,R]=strtok(R, '/');
% [bird,R]=strtok(R, '/');
% %[T,R]=strtok(R, '/');
% %[T,R]=strtok(R, '/');
% [site,R]=strtok(R,'/');
% 
%  [depth,R]=strtok(R, '/');
%  %[T,R]=strtok(R, '/');
%  
% %[W,R]=strtok(R,'/');
% %[D]=strtok(R,'/');

outfname = sprintf('%s_%s_%s_%s_ncmionto',subjectID,pen,site,sort); %use for cells labeled by location


 
 
  eval (['save ' outfname ' basenov basetrain l_songperinc u_songperinc baselifesparse druglifesparse washlifesparse basebias drugbias washbias']);
 
 !cp *ss*_ncmionto.mat ~/data/iontoncm/mean/.

 
 


%end



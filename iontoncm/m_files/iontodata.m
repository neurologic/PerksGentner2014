function []=iontodata(blocksize,blocktypes,sitetype)
% splits ionto data into separate files for each drug condition for further
% analysis
% blocksize = # reps in each blcok
% blocktypes = explains what drug was used in each block, in order 
%       1 = baseline
%       2 = gabazine
%       3 = gaba
%       4 = saclofen
%       5 = recovery
%       24 = gabazine + saclofen
%       23 = gabazine + gaba
% sitetype = what type of recording 
%       su = single unit
%       mu = multi unit
% 

%blocksize=10;
%blocktypes=[1 2 5];
%sitetype='mu';

doplot=0;
threshfactor=2;
plotlin=0;
binsize=50; %in msec
bps = 1000/binsize;
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};

cA4 = [0 .622 1.531 2.401 3.300  4.198  4.952  5.864  6.704 7.394  8.102  8.801 9.215];
cA2 = [0  .78  1.59  2.45  3.25  4.06  4.81  5.58  6.32 7.01  7.69  8.50 9.0756];
cB1 = [0 .731  1.698  2.707  3.456  4.228  5.247  5.996  6.536  7.030  7.614  8.300 8.916 9.571]; 
cB2 = [0  0.906  1.594  2.449  3.220  4.107  4.827  5.539  6.296  6.884  7.437  8.010 8.678];
cD1 = [0 0.844 1.749 2.921 4.312 5.578 5.971 6.710 7.423 8.185 8.481 9.140 9.491 10.122];
cD2 = [0 0.787 1.333 2.081 2.907 3.838 4.712 5.65 6.305 6.80 7.558 8.142 8.681 9.475 10.151];
cE1 = [0 0.637 1.440 2.193 2.739 3.605 4.042 4.543 5.087 5.871 6.547 7.443 8.208 8.903 9.350 9.950];
cE2 = [0 1.138 2.97 3.440 4.177 4.904 5.694 6.654 7.347 7.943 8.570 9.320 10.150];
scA4 = [0 .622 1.531 2.401 3.300  4.198  5.000];
scA2 = [0  .78  1.59  2.45  3.25  4.06  4.81  5.000];
scB1 = [0 .731  1.698  2.707  3.456  4.228  5.000]; 
scB2 = [0  0.906  1.594  2.449  3.220  4.107  4.827  5.000];
scD1 = [0 0.844 1.749 2.921 4.312 5.000];
scD2 = [0 0.787 1.333 2.081 2.907 3.838 4.712 5.000];
scE1 = [0 0.637 1.440 2.193 2.739 3.605 4.042 4.543 5.000];
scE2 = [0 1.138 2.97 3.440 4.177 5.000];

bird = {'st444'; 'st445'; 'st493'; 'st484'; 'st504'; 'st450'; 'st490'; 'st466'; 'st536'; 'st590'};
trainset = {'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'; 'E1' 'E2' 'A2' 'A4'; 'E1' 'E2' 'A2' 'A4'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'};
goset = {'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'; 'E1' 'E2'; 'E1' 'E2'; 'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'; 'A2' 'A4'; 'D1' 'D2'};
nogoset = { 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'; 'A2' 'A4'; 'A2' 'A4'; 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'; 'B1' 'B2'; 'E1' 'E2'};
novset = { 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'; 'B1' 'B2' 'D1' 'D2'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'};

%get trials for baseline condition
baserepsmax=find(blocktypes==1)*blocksize;
baserepsmin=baserepsmax-(blocksize-1);
basereps=baserepsmin:baserepsmax;

%get trials for gabazine condition
gzrepsmax=find(blocktypes==2)*blocksize;
gzrepsmin=gzrepsmax-(blocksize-1);
gzreps=gzrepsmin:gzrepsmax;

%get trials for recovery condition
recrepsmax=find(blocktypes==5)*blocksize;
recrepsmin=recrepsmax-(blocksize-1);
recreps=recrepsmin:recrepsmax;

%get trials for gaba condition if they exist
if isempty(find(blocktypes==3))==0
    gbrepsmax=find(blocktypes==3)*blocksize;
    gbrepsmin=gbrepsmax-(blocksize-1);
    gbreps=gbrepsmin:gbrepsmax;
end
%get trials for saclofen condition if they exist
if isempty(find(blocktypes==4))==0
    sacrepsmax=find(blocktypes==4)*blocksize;
    sacrepsmin=sacrepsmax-(blocksize-1);
    sacreps=sacrepsmin:sacrepsmax;
end
% get trials for gabazine+saclofen condition if exists 
if isempty(find(blocktypes==24))==0
    gzsacrepsmax=find(blocktypes==24)*blocksize;
    gzsacrepsmin=gzsacrepsmax-(blocksize-1);
    gzsacreps=gzsacrepsmin:gzsacrepsmax;
    
end
% get trials for gabazine+gaba condition if exists 
if isempty(find(blocktypes==23))==0
    gzgbrepsmax=find(blocktypes==23)*blocksize;
    gzgbrepsmin=gzgbrepsmax-(blocksize-1);
    gzgbreps=gzgbrepsmin:gzgbrepsmax;
    
end

pretime=2; %in sec
posttime=2; % in sec
files = dir('*dbmean_toe.txt');
allspikes=[];
trainstims=[];
novstims=[];
gostims=[];
nogostims=[];
Psth.Base=[];
Psth.Gz=[];
Psth.Rec=[];

x=1;
for i=1:length(files)
    toelist = files(i).name;
    %toelist = char(files(i));
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);  
    subindx=strcmpi(bird,subjectID);
    [stim,yt] = strtok(stimfile, '_');   
    stim = strtok(stim, '.');
    stim2=stim;
    if isempty(strfind(yt,'short'))==0
        stim2=strcat(stim2,'short');
    end
        
    stims{i}=stim2;
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
        
    stimdur = stimtimes(strmatch(stim2, stimnames, 'exact'));
    stimdurs(i)=stimdur;
    
    % get firing rates to individual motifs
    if stimdur<6
       cuts=eval(['sc' stim]);
       stimlength='s';
    else
       cuts = eval(['c' stim]); 
       stimlength='l'; 
    end
    
   

   % motifid = eval(['m' stim]);
    nummots(i)=size(cuts,2)-1;
    for y=1:size(cuts,2)-1;
            for l=1:nreps
                notoes = toes{l}{1}(find(toes{l}{:}>=cuts(y)));%only use spikes after stim starts
                notoes = notoes(find(notoes<=cuts(y+1))); %get rid of spikes after stim over
                motrate{i}(l,y)=length(notoes) / ( cuts(y+1)- cuts(y) );
                nummottoes(l,x)=length(notoes);
                motlength(l,x)=cuts(y+1)-cuts(y);
            end
            x=x+1;
    end
    
    xbins = linspace(0, stimdur, (stimdur)*bps);
    sponbins=linspace(-2,0,2*bps);
    allxbins{i}=xbins;
    numxbins(i)=length(xbins);
    psth=[];
    basesponpsth=[];
    gzsponpsth=[];
    recsponpsth=[];
    %othertoes = alltoes(find(alltoes>=0));
    %othertoes = othertoes(find(othertoes<=stimdur));
    for moo=1:nreps
            
            newtoes = toes{moo}{1}(find(toes{moo}{:}>=0));%only use spikes after stim starts
            newtoes = newtoes(find(newtoes<=stimdur)); %get rid of spikes after stim over 
            drivrate(moo,i) = length(newtoes)/(stimdur);
            spontoes=toes{moo}{1}(find(toes{moo}{:}<0));
            sponrate(moo,i) = length(spontoes )/2;
            
            if isempty(newtoes)==1
                tpsth=zeros(1,length(xbins));
            else 
                tpsth=histc(newtoes, xbins);
            end
            psth(moo,:) = tpsth';
            if isempty(spontoes)==1
                spsth=zeros(1,length(sponbins));
            else
                spsth=histc(spontoes,sponbins);
            end
            spsth=spsth./(binsize/1000);
            sponpsth(moo)=mean(spsth);
            
            bsponpsth(moo,:)=spsth';
            
            isis = diff(newtoes);
            isis = isis(find(isis>0));
            allisis{moo,i}=isis;
            numtoes(moo,i)=length(newtoes);
            
            spsth=reshape(spsth,1,length(sponbins));
            if isempty(find(basereps==moo))==0
               basesponpsth=[basesponpsth; spsth]; 
            end
            if isempty(find(gzreps==moo))==0
               gzsponpsth=[gzsponpsth; spsth]; 
            end
            if isempty(find(recreps==moo))==0
               recsponpsth=[recsponpsth; spsth]; 
            end
    end
    
        Psth.Base = [Psth.Base mean(psth(basereps,:))];
        Psth.Gz = [Psth.Gz mean(psth(gzreps,:)) ];
        Psth.Rec = [Psth.Rec mean(psth(recreps,:))];
    
    allpsth{i}=psth;
    allsponpsth(:,i)=sponpsth';
    
    allbsponpsth{i}=bsponpsth;
    BaseSponPsth(i,:)=mean(basesponpsth);
    GzSponPsth(i,:)=mean(gzsponpsth);
    RecSponPsth(i,:)=mean(recsponpsth);
end

[RBaseGzCorr,p]=corrcoef(Psth.Base, Psth.Gz);

basemotrate=[];
gzmotrate=[];
recmotrate=[];
tbasemot=[];
tgzmot=[];
trecmot=[];
nbasemot=[];
ngzmot=[];
nrecmot=[];
basemotreprate=[];
gzmotreprate=[];
recmotreprate=[];
for n=1:length(motrate)
   basemotrate=[basemotrate mean(motrate{n}(basereps,:))]; 
   gzmotrate=[gzmotrate mean(motrate{n}(gzreps,:))]; 
   recmotrate=[recmotrate mean(motrate{n}(recreps,:))]; 
   
   basemotreprate=[basemotreprate (motrate{n}(basereps,:))]; 
   gzmotreprate=[gzmotreprate (motrate{n}(gzreps,:))]; 
   recmotreprate=[recmotreprate (motrate{n}(recreps,:))]; 
   
    
   if isempty(find(trainstims==n))==0
        tbasemot=[tbasemot mean(motrate{n}(basereps,:))];
        tgzmot=[tgzmot mean(motrate{n}(gzreps,:))];
        trecmot=[trecmot mean(motrate{n}(recreps,:))];
   else
       nbasemot=[nbasemot mean(motrate{n}(basereps,:))];
        ngzmot=[ngzmot mean(motrate{n}(gzreps,:))];
        nrecmot=[nrecmot mean(motrate{n}(recreps,:))];        
   end
end

tgzpermot=(100*(tgzmot-tbasemot)./tbasemot);
ngzpermot=(100*(ngzmot-nbasemot)./nbasemot);


sponrates=mean(sponrate);
avgrate=mean(drivrate);
muspon = mean(sponrates);
threshrate=mean(sponrates)+std(sponrates);
numsigstims=length(find(avgrate>threshrate));
psigstims=numsigstims./length(stims);

depth=str2double(strtok(site,'Z'));
if mean(stimdurs) > 7
   stimlength='l'; 
else 
    stimlength='s';
end


basemotlength=motlength(basereps,:);
gzmotlength=motlength(gzreps,:);
recmotlength=motlength(recreps,:);

basesponrate=mean(sponrate(basereps,:));
baserate=mean(drivrate(basereps,:));
basereprate=drivrate(basereps,:);
basenumtoes=numtoes(basereps,:);
basenummottoes=nummottoes(basereps,:);
%if isempty(find(baserate==0))==0;
%    baserate(find(baserate==0))=0.01;
%end
basethresh=mean(basesponrate)+threshfactor*std(basesponrate);
basesuppthresh=mean(basesponrate)-threshfactor*std(basesponrate);
basesigstims=length(find(baserate>basethresh));
psigbasestims=basesigstims./length(stims);
zbaserate=zscore(baserate);
zbasetrain=zbaserate(trainstims);
zbasenov=zbaserate(novstims);
basenum=(sum(baserate./(length(baserate)))).^2;
baseden=sum((baserate.^2)./(length(baserate)));
baselifesparse=1-(basenum/baseden);
basetrain=mean(baserate(trainstims));
basenov=mean(baserate(novstims));
basebias=(basenov-basetrain)./(basenov+basetrain);
TBasePSigMots=length(find(tbasemot>basethresh))./length(tbasemot);
NBasePSigMots=length(find(nbasemot>basethresh))./length(nbasemot);
TBasePSuppMots=length(find(tbasemot<basesuppthresh))./length(tbasemot); 
NBasePSuppMots=length(find(nbasemot<basesuppthresh))./length(nbasemot); 
baseMotLifeSparse=LifetimeSparse(basemotrate);
baseTrainMotLifeSparse=LifetimeSparse(tbasemot);
baseNovelMotLifeSparse=LifetimeSparse(nbasemot);

sbasebinfr=mean(allsponpsth(basereps,:));
%binthresh=mean(mean(BaseSponPsth))+10*std(mean(BaseSponPsth));
binthresh=2*max(max(allsponpsth(basereps,:)));
gzbinthresh=2*max(max(allsponpsth(gzreps,:)));
%gzbinthresh=mean(mean(GzSponPsth))+10*std(mean(GzSponPsth));
allbasebinfr=[];
allgzbinfr=[];
binperinc=[];
AllBaseTrainBinFr=[];
AllBaseNovBinFr=[];
AllGzTrainBinFr=[];
AllGzNovBinFr=[];
for j=1:length(stims)
   currpsth=allpsth{j};

   basebinfr=mean(currpsth(basereps,:)./(binsize/1000));
   lbasebin(j)=length(basebinfr);
   gzbinfr=mean(currpsth(gzreps,:)./(binsize/1000));
   allbasebinfr=[allbasebinfr basebinfr];
   allgzbinfr=[allgzbinfr gzbinfr];
   if isempty(find(j==trainstims))==1
      AllBaseTrainBinFr=[AllBaseTrainBinFr basebinfr]; 
      AllGzTrainBinFr=[AllGzTrainBinFr gzbinfr];
   else
      AllBaseNovBinFr=[AllBaseNovBinFr basebinfr]; 
      AllGzNovBinFr=[AllGzNovBinFr gzbinfr];
   end
   basenoresp=find(basebinfr<=binthresh);
   nbasenoresp(j)=length(basenoresp);
   tgzbinfr=gzbinfr(basenoresp);
    gzuncov{j}=tgzbinfr(find(tgzbinfr>gzbinthresh));
    ngzuncov(j)=length(find(tgzbinfr>gzbinthresh));
    if isempty(find(tgzbinfr>gzbinthresh))==1
        NGzPeaksUncov(j)=0;
    else
        NGzPeaksUncov(j)=length(find(diff(find(tgzbinfr>gzbinthresh))>1))+1;
    end %gzbinfr(find(gzbinfr==0))=1;
  % basebinfr(find(basebinfr==0))=1;
  % binperinc=[binperinc 100*((gzbinfr-basebinfr)./basebinfr)];
   basebins=linspace(min(basebinfr),max(basebinfr),20);
   [counts,bin]=histc(basebinfr,basebins);
    for k=1:length(counts)
       binind=find(bin==k);
       gztempfr=gzbinfr(binind);
       [gzcounts]=histc(gztempfr,basebins);
       if sum(gzcounts>0)
       gzcounts=gzcounts./sum(gzcounts);
       end
       gzmat(k,:)=gzcounts;
    end
   % if plotlin%
%       figure;
%        imagesc(basebins,basebins,flipud(gzmat));
%        yvals=get(gca,'ytick');
%        set(gca,'yticklabel',fliplr(yvals))
%        xlabel('Gabazine Firing Rate')
%         ylabel('Baseline Firing Rate')
%        colormap(hot)


 %   end
    allgzmat{j}=gzmat;
end

pbasenoresp=sum(nbasenoresp)/sum(lbasebin);
pgzuncov=sum(ngzuncov)/sum(nbasenoresp);

if plotlin
   plot(allbasebinfr,allgzbinfr,'ok')
   xlabel('Binned Baseline Firing Rate')
   ylabel('Binned Gabazine Firing Rate')
    
    
   
   
   
end
baseisis=[];
gzisis=[];
recisis=[];
for l=1:blocksize
   for m=1:length(stims)
        baseisis=[baseisis; allisis{basereps(l),m} ];
        gzisis=[gzisis; allisis{gzreps(l),m} ];
        recisis=[recisis; allisis{recreps(l),m} ];
   end
end
baseisi=mean(baseisis);
gzisi=mean(gzisis);
recisi=mean(recisis);

  if doplot
     figure
     hist(allgzbinfr-allbasebinfr,50);
     xlabel('Magnitude of Firing Rate Increase with Gabazine (spikes/sec)')
     ylabel(' # 50ms bins')
     h = findobj(gca,'Type','patch');
     set(h,'FaceColor','k')
      
     figure
     tallgzbinfr=allgzbinfr;
     tallbasebinfr=allbasebinfr;
     tallgzbinfr(find(tallgzbinfr==0))=1;
     tallbasebinfr(find(tallbasebinfr==0))=1;
     hist(100*(tallgzbinfr-tallbasebinfr)./tallbasebinfr,50);
     xlabel('% Increase in Firing Rate  with Gabazine (spikes/sec)')
     ylabel(' # 50ms bins')
     h = findobj(gca,'Type','patch');
     set(h,'FaceColor','k')
      
      
      
  end
        
gzsponrate=mean(sponrate(gzreps,:));
gzsponper=100*((gzsponrate-basesponrate)./(basesponrate));
gzrate=mean(drivrate(gzreps,:));
gzreprate=drivrate(gzreps,:);
gznumtoes=numtoes(gzreps,:);
gznummottoes=nummottoes(gzreps,:);
gzthresh=mean(gzsponrate)+threshfactor*std(gzsponrate);
gzsuppthresh=mean(gzsponrate)-threshfactor*std(gzsponrate);
TGzPSigMots=length(find(tgzmot>gzthresh))./length(tgzmot);
NGzPSigMots=length(find(ngzmot>gzthresh))./length(ngzmot);

TGzPSuppMots=length(find(tgzmot<gzsuppthresh))./length(tgzmot); 
NGzPSuppMots=length(find(ngzmot<gzsuppthresh))./length(ngzmot); 
gzsigstims=length(find(gzrate>gzthresh));
psiggzstims=gzsigstims./length(stims);
%gzper=zscore((gzrate-mean(baserate))./mean(baserate));
Tempbaserate=baserate;
Tempbaserate(find(Tempbaserate==0))=0.01;
Tempgzrate=gzrate;
Tempgzrate(find(Tempgzrate==0))=0.01;
gzper=100*(gzrate-Tempbaserate)./(Tempbaserate);
traingzper=gzper(trainstims);
novgzper=gzper(novstims);
gzdiff=gzrate-baserate;
traingzdiff=gzdiff(trainstims);
novgzdiff=gzdiff(novstims);

gogzper=gzper(gostims);
nogogzper=gzper(nogostims);
[hyp,gzperpval]=ttest2(traingzper,novgzper);
%traingzper=100*((mean(gzrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
%novgzper=100*((mean(gzrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
zgzrate=zscore(gzrate);
zgztrain=zgzrate(trainstims);
zgznov=zgzrate(novstims);
gznum=(sum(gzrate./(length(gzrate)))).^2;
gzden=sum((gzrate.^2)./(length(gzrate)));
gzlifesparse=1-(gznum/gzden);
gztrain=mean(gzrate(trainstims));
gznov=mean(gzrate(novstims));
gzbias=(gznov-gztrain)./(gznov+gztrain);
gzMotLifeSparse=LifetimeSparse(gzmotrate);
gzTrainMotLifeSparse=LifetimeSparse(tgzmot);
gzNovelMotLifeSparse=LifetimeSparse(ngzmot);

DiffMotLifeSparse = LifetimeSparse(gzmotrate-basemotrate);
DiffNovMotLifeSparse = LifetimeSparse(ngzmot-nbasemot);
DiffTrainMotLifeSparse = LifetimeSparse(tgzmot-tbasemot);

AllRate = [baserate gzrate];
zAllRate = zscore(AllRate);
AllBaseTrainStims=[trainstims];
AllGzTrainStims=[length(stims)+trainstims];
AllBaseNovStims=[novstims ];
AllGzNovStims=[length(stims)+novstims];


recsponrate=mean(sponrate(recreps,:));
recreprate=drivrate(recreps,:);
recnumtoes=numtoes(recreps,:);
recnummottoes=nummottoes(recreps,:);
recsponper=100*((recsponrate-basesponrate)./(basesponrate));
recrate=mean(drivrate(recreps,:));
recthresh=mean(recsponrate)+threshfactor*std(recsponrate);
recsuppthresh=mean(recsponrate)-threshfactor*std(recsponrate);
TRecPSigMots=length(find(trecmot>recthresh))./length(trecmot);
NRecPSigMots=length(find(nrecmot>recthresh))./length(nrecmot);
TRecPSuppMots=length(find(trecmot<recsuppthresh))./length(trecmot); 
NRecPSuppMots=length(find(nrecmot<recsuppthresh))./length(nrecmot); 


AllRatewRec = [baserate gzrate recrate];
zAllRatewRec = zscore(AllRatewRec);
AllwRecBaseTrainStims = trainstims;
AllwRecGzTrainStims = length(stims)+trainstims;
AllwRecBaseNovStims = novstims;
AllwRecGzNovStims = length(stims)+novstims;
AllwRecRecTrainStims = trainstims+length(stims)+length(stims);
AllwRecRecNovStims = novstims+length(stims)+length(stims);


recsigstims=length(find(recrate>recthresh));
psigrecstims=recsigstims./length(stims);
recper=100*(recrate-mean(baserate))./mean(baserate);
trainrecper=recper(trainstims);
novrecper=recper(novstims);
%trainrecper=100*((mean(recrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
%novrecper=100*((mean(recrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
zrecrate=zscore(recrate);
zrectrain=zrecrate(trainstims);
zrecnov=zrecrate(novstims);
recnum=(sum(recrate./(length(recrate)))).^2;
recden=sum((recrate.^2)./(length(recrate)));
reclifesparse=1-(recnum/recden);
rectrain=mean(recrate(trainstims));
recnov=mean(recrate(novstims));
recbias=(recnov-rectrain)./(recnov+rectrain);
recMotLifeSparse=LifetimeSparse(recmotrate);
recTrainMotLifeSparse=LifetimeSparse(trecmot);
recNovelMotLifeSparse=LifetimeSparse(nrecmot);



R = 1000;                % number of simulations

% ask whether Baseline Bias value greater than expected by chance
for r = 1:R
%shuffle real rates on each trialas learned and unfam
    permnorm=[];
    for n=1:length(basereps)
        temprate=basereprate(n,:);
        temprate=temprate(randperm(length(temprate)));
        permnorm=[permnorm; temprate];

    end
    mpermnorm=mean(permnorm);
    permnovrate2 = mean(mpermnorm(novstims));
    permtrainrate2 = mean(mpermnorm(trainstims));
    scoreR = (permnovrate2-permtrainrate2)/(permnovrate2+permtrainrate2);
    SMATRat4(r,:) = scoreR; 
end

UCIBaseRate = prctile(SMATRat4(:), 95);

if basebias >=0;
    PCNTILBaseRate = 1 - (length(find(SMATRat4(:)<basebias))/length(SMATRat4(:))) ;
else
    PCNTILBaseRate = 1 - (length(find(SMATRat4(:)>basebias))/length(SMATRat4(:))) ;
end

% ask whether Gabazine Bias value greater than expected by chance
for r = 1:R
%shuffle real rates on each trialas learned and unfam
    permnorm=[];
    for n=1:length(basereps)
        temprate=gzreprate(n,:);
        temprate=temprate(randperm(length(temprate)));
        permnorm=[permnorm; temprate];

    end
    mpermnorm=mean(permnorm);
    permnovrate2 = mean(mpermnorm(novstims));
    permtrainrate2 = mean(mpermnorm(trainstims));
    scoreR = (permnovrate2-permtrainrate2)/(permnovrate2+permtrainrate2);
    SMATRat4(r,:) = scoreR; 
end

UCIGzRate = prctile(SMATRat4(:), 95);

if gzbias >=0;
    PCNTILGzRate = 1 - (length(find(SMATRat4(:)<gzbias))/length(SMATRat4(:))) ;
else
    PCNTILGzRate = 1 - (length(find(SMATRat4(:)>gzbias))/length(SMATRat4(:))) ;
end

% Shuffle MEan FRs across stims
for r = 1:R
%shuffle mean firing rates among stimuli
permrate=baserate(randperm(length(baserate)));
permnovrate = mean(permrate(1:length(novstims)));
permtrainrate = mean(permrate(length(novstims)+1:length(stims)));
scoreR = (permnovrate-permtrainrate)/(permnovrate+permtrainrate);
SMATRat3(r,:) = scoreR; 



end

UCIRat3 = prctile(SMATRat3(:), 95);

if basebias >=0;
    PCNTILBaseRate2 = 1 - (length(find(SMATRat3(:)<basebias))/length(SMATRat3(:))) ;
else
    PCNTILBaseRate2 = 1 - (length(find(SMATRat3(:)>basebias))/length(SMATRat3(:))) ;
end

% Shuffle MEan FRs across stims
for r = 1:R
%shuffle mean firing rates among stimuli
permrate=gzrate(randperm(length(gzrate)));
permnovrate = mean(permrate(1:length(novstims)));
permtrainrate = mean(permrate(length(novstims)+1:length(stims)));
scoreR = (permnovrate-permtrainrate)/(permnovrate+permtrainrate);
SMATRat3(r,:) = scoreR; 



end

UCIRat3 = prctile(SMATRat3(:), 95);

if gzbias >=0;
    PCNTILGzRate2 = 1 - (length(find(SMATRat3(:)<gzbias))/length(SMATRat3(:))) ;
else
    PCNTILGzRate2 = 1 - (length(find(SMATRat3(:)>gzbias))/length(SMATRat3(:))) ;
end


if isempty(find(blocktypes==3))==0
    gbsponrate=mean(sponrate(gbreps,:));
    gbsponper=100*((gbsponrate-basesponrate)./(basesponrate));
    gbrate=mean(drivrate(gbreps,:));
    gbthresh=mean(gbsponrate)+threshfactor*std(gbsponrate);
    gbsigstims=length(find(gbrate>gbthresh));
    psiggbstims=gbsigstims./length(stims);
    gbper=100*(gbrate-mean(baserate))./mean(baserate);
    traingbper=gbper(trainstims);
    novgbper=gbper(novstims);
    %traingbper=100*((mean(gbrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
    %novgbper=100*((mean(gbrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
    zgbrate=zscore(gbrate);
    zgbtrain=zgbrate(trainstims);
    zgbnov=zgbrate(novstims);
    gbnum=(sum(gbrate./(length(gbrate)))).^2;
    gbden=sum((gbrate.^2)./(length(gbrate)));
    gblifesparse=1-(gbnum/gbden);
    gbtrain=mean(gbrate(trainstims));
    gbnov=mean(gbrate(novstims));
    gbbias=(gbnov-gbtrain)./(gbnov+gbtrain);
end
if isempty(find(blocktypes==4))==0
    sacsponrate=mean(sponrate(sacreps,:));
    sacsponper=100*((sacsponrate-basesponrate)./(basesponrate));
    sacrate=mean(drivrate(sacreps,:));
    sacthresh=mean(sacsponrate)+threshfactor*std(sacsponrate);
    sacsigstims=length(find(sacrate>sacthresh));
    psigsacstims=sacsigstims./length(stims);
    sacper=100*(sacrate-mean(baserate))./mean(baserate);
    trainsacper=sacper(trainstims);
    novsacper=sacper(novstims);
    %trainsacper=100*((mean(sacrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
    %novsacper=100*((mean(sacrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
    zsacrate=zscore(sacrate);
    zsactrain=zsacrate(trainstims);
    zsacnov=zsacrate(novstims);
    sacnum=(sum(sacrate./(length(sacrate)))).^2;
    sacden=sum((sacrate.^2)./(length(sacrate)));
    saclifesparse=1-(sacnum/sacden);
    sactrain=mean(sacrate(trainstims));
    sacnov=mean(sacrate(novstims));
    sacbias=(sacnov-sactrain)./(sacnov+sactrain);
end
if isempty(find(blocktypes==24))==0
    gzsacsponrate=mean(sponrate(gzsacreps,:));
    gzsacsponper=100*((gzsacsponrate-basesponrate)./(basesponrate));
    gzsacrate=mean(drivrate(gzsacreps,:));
    gzsacthresh=mean(gzsacsponrate)+threshfactor*std(gzsacsponrate);
    gzsacsigstims=length(find(gzsacrate>gzsacthresh));
    psiggzsacstims=gzsacsigstims./length(stims);
    gzsacper=100*(gzsacrate-mean(baserate))./mean(baserate);
    traingzsacper=gzsacper(trainstims);
    novgzsacper=gzsacper(novstims);
    traingzsacper=100*((mean(gzsacrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
    novgzsacper=100*((mean(gzsacrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
    zgzsacrate=zscore(gzsacrate);
    zgzsactrain=zgzsacrate(trainstims);
    zgzsacnov=zgzsacrate(novstims);
    gzsacnum=(sum(gzsacrate./(length(gzsacrate)))).^2;
    gzsacden=sum((gzsacrate.^2)./(length(gzsacrate)));
    gzsaclifesparse=1-(gzsacnum/gzsacden);
    gzsactrain=mean(gzsacrate(trainstims));
    gzsacnov=mean(gzsacrate(novstims));
    gzsacbias=(gzsacnov-gzsactrain)./(gzsacnov+gzsactrain);
end
if isempty(find(blocktypes==23))==0
    gzgbsponrate=mean(sponrate(gzgbreps,:));
    gzgbsponper=100*((gzgbsponrate-basesponrate)./(basesponrate));
    gzgbrate=mean(drivrate(gzgbreps,:));
    gzgbthresh=mean(gzgbsponrate)+threshfactor*std(gzgbsponrate);
    gzgbsigstims=length(find(gzgbrate>gzgbthresh));
    psiggzgbstims=gzgbsigstims./length(stims);
    gzgbper=100*(gzgbrate-mean(baserate))./mean(baserate);
    traingzgbper=gzgbper(trainstims);
    novgzgbper=gzgbper(novstims);
    %traingzgbper=100*((mean(gzgbrate(trainstims))-mean(baserate(trainstims)))/(mean(baserate(trainstims))));
    %novgzgbper=100*((mean(gzgbrate(novstims))-mean(baserate(novstims)))/(mean(baserate(novstims))));
    zgzgbrate=zscore(gzgbrate);
    zgzgbtrain=zgzgbrate(trainstims);
    zgzgbnov=zgzgbrate(novstims);
    gzgbnum=(sum(gzgbrate./(length(gzgbrate)))).^2;
    gzgbden=sum((gzgbrate.^2)./(length(gzgbrate)));
    gzgblifesparse=1-(gzgbnum/gzgbden);
    gzgbtrain=mean(gzgbrate(trainstims));
    gzgbnov=mean(gzgbrate(novstims));
    gzgbbias=(gzgbnov-gzgbtrain)./(gzgbnov+gzgbtrain);
end

clear files i 




outfname = sprintf('%s_%s_%s_%s_iontodata',subjectID,pen,site,sort); %use for cells labeled by location
save(outfname);
 
 !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/all/.
if strcmp(sitetype,'mu')==1
 !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/mu/.
else
 !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/su/.   
end

 
    

    
end

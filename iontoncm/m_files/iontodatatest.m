%function []=iontodata(blocksize,blocktypes,sitetype)
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

blocksize=10;
blocktypes=[1 2 5 ];
%sitetype='mu';
plotlin=0;
binsize=50; %in msec
bps = 1000/binsize;
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};


bird = {'st444'; 'st445'; 'st493'; 'st484'; 'st504'; 'st450'; 'st490'; 'st466'};
trainset = {'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'; 'E1' 'E2' 'A2' 'A4'; 'E1' 'E2' 'A2' 'A4'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'};
goset = {'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'; 'E1' 'E2'; 'E1' 'E2'; 'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'};
nogoset = { 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'; 'A2' 'A4'; 'A2' 'A4'; 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'};
novset = { 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'; 'B1' 'B2' 'D1' 'D2'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'};

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
    
    xbins = linspace(0, stimdur, (stimdur)*bps);
    sponbins=linspace(-2,0,2*bps);
    allxbins{i}=xbins;
    numxbins(i)=length(xbins);
    psth=[];
    %othertoes = alltoes(find(alltoes>=0));
    %othertoes = othertoes(find(othertoes<=stimdur));
    for moo=1:nreps
            
            newtoes = toes{moo}{1}(find(toes{moo}{:}>=0));%only use spikes after stim starts
            newtoes = newtoes(find(newtoes<=stimdur)); %get rid of spikes after stim over 
            drivrate(moo,i) = length(newtoes)/(stimdur);
            spontoes=toes{moo}{1}(find(toes{moo}{:}<0));
            sponrate(moo,i) = length(spontoes )/2;
            
            spikes{i,moo}=newtoes;    
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
           
            
    end
    allpsth{i}=psth;
    allsponpsth(:,i)=sponpsth';
    
end
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

threshfactor=1; %in SDs of spontaneous rate
basesponrate=mean(sponrate(basereps,:));
baserate=mean(drivrate(basereps,:));
%if isempty(find(baserate==0))==0;
%    baserate(find(baserate==0))=0.01;
%end
basethresh=mean(basesponrate)+threshfactor*std(basesponrate);
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

sbasebinfr=mean(allsponpsth(basereps,:));
binthresh=mean(sbasebinfr)+std(sbasebinfr);

binperinc=[];
for j=1:length(stims)
   currpsth=allpsth{j};

   basebinfr=mean(currpsth(basereps,:)./(binsize/1000));
   gzbinfr=mean(currpsth(gzreps,:)./(binsize/1000));
   basenoresp=find(basebinfr<binthresh);
   tgzbinfr=gzbinfr(basenoresp);
    gzuncov{j}=tgzbinfr(find(tgzbinfr>binthresh));
   %gzbinfr(find(gzbinfr==0))=1;
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
    if plotlin
       figure;
       imagesc(basebins,basebins,flipud(gzmat));
       yvals=get(gca,'ytick');
       set(gca,'yticklabel',fliplr(yvals))
       xlabel('Gabazine Firing Rate')
        ylabel('Baseline Firing Rate')
       colormap(hot)
    end
    allgzmat{j}=gzmat;
end



gzsponrate=mean(sponrate(gzreps,:));
gzsponper=100*((gzsponrate-basesponrate)./(basesponrate));
gzrate=mean(drivrate(gzreps,:));
gzthresh=mean(gzsponrate)+threshfactor*std(gzsponrate);
gzsigstims=length(find(gzrate>gzthresh));
psiggzstims=gzsigstims./length(stims);
%gzper=zscore((gzrate-mean(baserate))./mean(baserate));
gzper=100*(gzrate-mean(baserate))./mean(baserate);
traingzper=gzper(trainstims);
novgzper=gzper(novstims);
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

recsponrate=mean(sponrate(recreps,:));
recsponper=100*((recsponrate-basesponrate)./(basesponrate));
recrate=mean(drivrate(recreps,:));
recthresh=mean(recsponrate)+threshfactor*std(recsponrate);
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
 
%  !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/all/.
% if strcmp(sitetype,'mu')==1
%  !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/mu/.
% else
%  !cp *ss*_iontodata.mat ~/data/iontoncm/splitbydrug/su/.   
% end
% 
%  
%     
% 
%     
%end

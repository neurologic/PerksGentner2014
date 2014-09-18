%function []=iontodatamean(sitetype)
% sitetype = what type of recordings to include in analysis 
%       su = single unit
%       mu = multi unit
%       both = both su and mu
%%
clear all
close all

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

sitetype='both';
doplot=0;

if strcmp(sitetype,'both')==1
    
    cd ../splitbydrug/all/ 
elseif strcmp(sitetype,'su')==1
    cd ../splitbydrug/su/
elseif strcmp(sitetype,'mu')==1
    cd ../splitbydrug/mu/
end
gbindex=[];
mgbper=[];
mgbsponper=[];
mtraingbper=[];
mnovgbper=[];
mgblsparse=[]; 
mzgbtrain=[];
mzgbnov=[];
mpsiggb=[];
mgbbias=[];
mgbrate=[];
sacindex=[];
msacper=[];
msacsponper=[];
mtrainsacper=[];
mnovsacper=[];
msaclsparse=[];
mzsactrain=[];
mzsacnov=[];
mpsigsac=[];
msacbias=[];
msacrate=[]; 
gzsacindex=[];
gzsacindex=[];
mgzsacper=[];
mgzsacsponper=[];
mtraingzsacper=[];
mnovgzsacper=[];
mgzsaclsparse=[]; 
mzgzsactrain=[];
mzgzsacnov=[];
mpsiggzsac=[];
mgzsacbias=[];
mgzsacrate=[];
        
gzgbindex=[];
gzgbindex=[];
mgzgbper=[];
mgzgbsponper=[];
mtraingzgbper=[];
mnovgzgbper=[];
mgzgblsparse=[]; 
mzgzgbtrain=[];
mzgzgbnov=[];
mpsiggzgb=[];
mgzgbbias=[];
mgzgbrate=[];
 
gzinc=[];
gzdec=[];
gbinc=[];
gbdec=[];
sacinc=[];
sacdec=[];
gzsacinc=[];
gzsacdec=[];
allbasemot=[];
allgzinc=[];

SigGzLearnPref=[];
SigGzUnfamPref=[];
SigBaseLearnPref=[];
SigBaseUnfamPref=[];

SigGzLearnMotPref=[];
SigGzUnfamMotPref=[];
SigBaseLearnMotPref=[];
SigBaseUnfamMotPref=[];

AllSitesGzIncBin=[];
AllSitesGzPerIncBin=[];
AllSiteBinBaseFR=[];
AllSiteBinGzFR=[];

%BinFRBins=linspace(-20,20,30);
BinFRBins=-100:50:500;
g=1;
f=1;
datafiles=dir('*iontodata.mat');
for i=1:length(datafiles)
   load(datafiles(i).name)
   birdID{i}=subjectID;
   sitetypes{i}=sitetype;
   allblocktypes{i}=blocktypes;
   
   RCorrPsth(i)=RBaseGzCorr(2);
   if stimlength=='l';
        LongMBaseMotRate{f}=basemotrate;
          LongMGzMotRate{f}=gzmotrate;
          f=f+1;
   elseif stimlength=='s';
        ShortMBaseMotRate{g}=basemotrate;
        ShortMGzMotRate{g}=gzmotrate;
   g=g+1;
   end
   stimlengths{i}=stimlength;
   depths(i)=depth;
   mbaselsparse(i)=(baselifesparse);
   basenum=(sum(basemotrate./(length(basemotrate)))).^2;
    baseden=sum((basemotrate.^2)./(length(basemotrate)));
    mbaselsparsemot(i)=1-(basenum/baseden);
   MbaseMotLifeSparse(i) = LifetimeSparsness(basemotrate);
   MbaseTrainMotLifeSparse(i) = baseTrainMotLifeSparse;
   MbaseNovelMotLifeSparse(i) = baseNovelMotLifeSparse;   
   mzbasetrain(i)=mean(zbasetrain);
   mzbasenov(i)=mean(zbasenov);
   MZBaseTrain2(i)=mean(zAllRate(AllBaseTrainStims));
   MZBaseNov2(i)=mean(zAllRate(AllBaseNovStims));
   MZGzTrain2(i)=mean(zAllRate(AllGzTrainStims));
   MZGzNov2(i)=mean(zAllRate(AllGzNovStims));
   
   MDiffMotLifeSparse(i)=DiffMotLifeSparse;
   MDiffNovMotLifeSparse(i)=DiffNovMotLifeSparse;
   MDiffTrainMotLifeSparse(i)=DiffTrainMotLifeSparse;
   
   %zscore for train and nov taken across base, gz and recov
   MZBaseTrain3(i)=mean(zAllRatewRec(AllwRecBaseTrainStims));
   MZBaseNov3(i)=mean(zAllRatewRec(AllwRecBaseNovStims));
   MZGzTrain3(i)=mean(zAllRatewRec(AllwRecGzTrainStims));   
   MZGzNov3(i)=mean(zAllRatewRec(AllwRecGzNovStims));      
   MZRecTrain3(i)=mean(zAllRatewRec(AllwRecRecTrainStims));      
   MZRecNov3(i)=mean(zAllRatewRec(AllwRecRecNovStims));      
   
   %Zscore taken across Base Gz and Rec
   zAllRate2=zscore([baserate gzrate recrate]);
   MZBaseRate(i)=mean(zAllRate2(1:8));
   MZGzRate(i)=mean(zAllRate2(9:16));
   MZRecRate(i)=mean(zAllRate2(17:24));
   
   %Zscore for Spon Rate taken across Base Gz Rec
   zAllSponRate2=zscore([basesponrate gzsponrate recsponrate]);
   MZBaseSponRate(i)=mean(zAllSponRate2(1:8));
   MZGzSponRate(i)=mean(zAllSponRate2(9:16));
   MZRecSponRate(i)=mean(zAllSponRate2(17:24));
       
   MTBasePMotExc(i) = mean(TBasePSigMots);
   MTGzPMotExc(i) = mean(TGzPSigMots);
   MTRecPMotExc(i) = mean(TRecPSigMots);
   MTBasePMotSup(i) = mean(TBasePSuppMots);
   MTGzPMotSup(i) = mean(TGzPSuppMots);
   MTRecPMotSup(i) = mean(TRecPSuppMots);
   MNBasePMotExc(i) = mean(NBasePSigMots);
   MNGzPMotExc(i) = mean(NGzPSigMots);
   MNRecPMotExc(i) = mean(NRecPSigMots);
   MNBasePMotSup(i) = mean(NBasePSuppMots);
   MNGzPMotSup(i) = mean(NGzPSuppMots);
   MNRecPMotSup(i) = mean(NRecPSuppMots);
   
   MeanPCNTILBaseRate(i)=PCNTILBaseRate;
   MeanPCNTILGzRate(i) = PCNTILGzRate;
   MeanPCNTILBaseRate2(i)=PCNTILBaseRate2;
   MeanPCNTILGzRate2(i) = PCNTILGzRate2;
   
   ZGzDiff=zscore(gzdiff);
   MZTrainDiff(i)=mean(ZGzDiff(trainstims));
   MZNovDiff(i)=mean(ZGzDiff(novstims));
   
   mpsignbase(i)=psigbasestims;
   mbasebias(i)=basebias;
   mbaserate(i)=mean(baserate);
   mbasesponrate(i)=mean(basesponrate);
   mbasetrain(i)=basetrain;
   mbasenov(i)=basenov;
   mgztrain(i)=gztrain;
   mgznov(i)=gznov;
   mbaseisi(i)=baseisi;
   allpbasenoresp(i)=pbasenoresp;
   mtbasemot(i)=mean(tbasemot);
   mnbasemot(i)=mean(nbasemot);
   tzeromots=intersect(find(tbasemot<basethresh),find(tgzmot<gzthresh));
   nzeromots=intersect(find(nbasemot<basethresh),find(ngzmot<gzthresh));
   tbasemot2=tbasemot;
   nbasemot2=nbasemot;
   tbasemot2(tzeromots)=[];
   nbasemot2(nzeromots)=[];
   mtbasemot2(i)=mean(tbasemot);% omit motifs with zero FR in both gz and base conditions
   mnbasemot2(i)=mean(nbasemot);
   if min(baserate)==0
       baserate(find(baserate==min(baserate)))=0.01;
   end
   mbasedyn(i)=max(baserate)./min(baserate);
   
   mgzper(i)=mean(gzper);
   mgzsponper(i)=mean(gzsponper);
   
   mtraingzper(i)=mean(traingzper);
   mnovgzper(i)=mean(novgzper);
   mgzdiff(i)=mean(gzdiff);
   GzZDiff=zscore(gzdiff);
   MTrainGzZDiff(i)=mean(GzZDiff(trainstims));
   MNovGzZDiff(i)=mean(GzZDiff(novstims));
   mtraingzdiff(i)=mean(traingzdiff);
   mnovgzdiff(i)=mean(novgzdiff);
   mgogzper(i)=mean(gogzper);
   mnogogzper(i)=mean(nogogzper);
   mgzlsparse(i)=(gzlifesparse); 
   mzgztrain(i)=mean(zgztrain);
   mzgznov(i)=mean(zgznov);
   mpsiggz(i)=psiggzstims;
   mgzbias(i)=gzbias;
   mgzperpval(i)=gzperpval;
   mgzrate(i)=mean(gzrate);
   mgzsponrate(i)=mean(gzsponrate);
   mgzisi(i)=gzisi;
   allpgzuncov(i)=pgzuncov;
   if min(gzrate)==0
       gzrate(find(gzrate==min(gzrate)))=0.01;
   end
   mgzdyn(i)=max(gzrate)./min(gzrate);
   mgzinc(i)=mean(allgzbinfr-allbasebinfr);
   AllSiteBinFRBase{i} = allbasebinfr;
   AllSiteBinFRGz{i} = allgzbinfr;
   
   AllSiteBinBaseFRCV(i) = std(allbasebinfr)/mean(allbasebinfr);
   AllSiteBinGzFRCV(i) = std(allgzbinfr)/mean(allgzbinfr);  
   AllSiteBinBaseTrainFRCV(i)=std(AllBaseTrainBinFr)/mean(AllBaseTrainBinFr);
   AllSiteBinBaseNovFRCV(i)=std(AllBaseNovBinFr)/mean(AllBaseNovBinFr);
   AllSiteBinGzTrainFRCV(i)=std(AllGzTrainBinFr)/mean(AllGzTrainBinFr);
   AllSiteBinGzNovFRCV(i)=std(AllGzNovBinFr)/mean(AllGzNovBinFr);
   
   AllSiteBinBaseFRSTD(i) = std(allbasebinfr);
   AllSiteBinGzFRSTD(i) = std(allgzbinfr); 
   AllSiteBinBaseFR = [AllSiteBinBaseFR allbasebinfr];
   AllSiteBinGzFR = [AllSiteBinGzFR allgzbinfr];
   AllSiteBinBaseTrainFRSTD(i)=std(AllBaseTrainBinFr);
   AllSiteBinBaseNovFRSTD(i)=std(AllBaseNovBinFr);
   AllSiteBinGzTrainFRSTD(i)=std(AllGzTrainBinFr);
   AllSiteBinGzNovFRSTD(i)=std(AllGzNovBinFr);
      
   allgzbinfr(find(allgzbinfr==0))=1;
   allbasebinfr(find(allbasebinfr==0))=1;
   AllSiteBinFRPerIncMat(i,:)=hist(100*((allgzbinfr-allbasebinfr)./allbasebinfr),BinFRBins);
   NumFRBins=length(allbasebinfr);
   BinFRNormMat(i,:)=ones(1,length(BinFRBins))+(NumFRBins-1);
   
   
   AllSitesGzIncBin = [AllSitesGzIncBin allgzbinfr-allbasebinfr];
   AllSitesGzPerIncBin = [AllSitesGzPerIncBin 100*(allgzbinfr-allbasebinfr)./allbasebinfr];
   AllSiteGzIncMax(i)=max(allgzbinfr-allbasebinfr);
   AllSiteGzIncMin(i)=min(allgzbinfr-allbasebinfr); 
   AllSiteGzPerIncMax(i)=max(100*(allgzbinfr-allbasebinfr)./allbasebinfr);
   AllSiteGzPerIncMin(i)=min(100*(allgzbinfr-allbasebinfr)./allbasebinfr); 
   medgzinc(i)=median(allgzbinfr-allbasebinfr);
   iqrgzinc(i)=iqr(allgzbinfr-allbasebinfr);
  
   
   
   mgzninc(i)=mean(100*(allgzbinfr-allbasebinfr)./allbasebinfr);
   medgzninc(i)=median(100*(allgzbinfr-allbasebinfr)./allbasebinfr);
   iqrgzninc(i)=iqr(100*(allgzbinfr-allbasebinfr)./allbasebinfr);
   allgzninc{i}=100*(allgzbinfr-allbasebinfr)./allbasebinfr;
   allgzincindex{i}=ones(1,length(allgzbinfr))+(i-1);
   [gzp,h]=signrank(baserate,gzrate);
   gzpvals(i)=gzp;
   if gzp<0.05
       if mean(gzrate)>mean(baserate)
           gzinc=[gzinc i];
       end
       if mean(gzrate)<mean(baserate)
           gzdec=[gzdec i];
       end
   end
   gztrainuncov=[];
   gznovuncov=[];
   for m=1:length(trainstims)
        gztrainuncov=[gztrainuncov gzuncov{trainstims(m)}];
        gznovuncov=[gznovuncov gzuncov{novstims(m)}];
   end   
   
   [h,p1]=ttest2(baserate(trainstims),baserate(novstims));
   PValBaseFR(i) = p1;
   if p1 < 0.05
       if mean(baserate(trainstims))> mean(baserate(novstims))
           SigBaseLearnPref = [SigBaseLearnPref i];
       else
           SigBaseUnfamPref = [SigBaseUnfamPref i];
       end
   end   
   [h,p2]=ttest2(gzrate(trainstims),gzrate(novstims));
   PValGzFR(i) = p2;
   if p2 < 0.05
       if mean(gzrate(trainstims))> mean(gzrate(novstims))
           SigGzLearnPref = [SigGzLearnPref i];
       else
           SigGzUnfamPref = [SigGzUnfamPref i];
       end
   end   
   
   [h,p3]=ttest2(tbasemot,nbasemot);
   PValBaseMotFR(i) = p3;
   if p3 < 0.05
       if mean(tbasemot)> mean(nbasemot)
           SigBaseLearnMotPref = [SigBaseLearnMotPref i];
       else
           SigBaseUnfamMotPref = [SigBaseUnfamMotPref i];
       end
   end   
   [h,p4]=ttest2(tgzmot,ngzmot);
   PValGzMotFR(i) = p4;
   if p4 < 0.05
       if mean(tgzmot)> mean(ngzmot)
           SigGzLearnMotPref = [SigGzLearnPref i];
       else
           SigGzUnfamMotPref = [SigGzUnfamPref i];
       end
   end   
   mpgzuncov(i)=pgzuncov;
   mpbasenoresp(i)=pbasenoresp;
   songnoresp=nbasenoresp./lbasebin;
   songnorespindex=find(songnoresp==1); %think about another measure of no resp e.g. >99%
   mnsongnoresp(i)=length(find(ngzuncov(songnorespindex)>0));
   % ngztrainuncov(i)=length(gztrainuncov);
  % ngznonvuncov(i)=length(gznovuncov);
   ngztrainuncov(i)=length(gztrainuncov)./sum(nbasenoresp(trainstims));
   ngznovuncov(i)=length(gznovuncov)./sum(nbasenoresp(novstims));
   NGzTrainPeaksUncov(i) = sum(NGzPeaksUncov(trainstims));
   NGzNovPeaksUncov(i) = sum(NGzPeaksUncov(novstims));
   NGzPeaks(i) = sum(NGzPeaksUncov);
   FGzPeaks(i) = sum(NGzPeaksUncov)./sum(stimdurs);
   FGzTrainPeaks(i) = sum(NGzPeaksUncov(trainstims))./sum(stimdurs(trainstims));
   FGzNovPeaks(i) = sum(NGzPeaksUncov(novstims))./sum(stimdurs(novstims));
   PGzTrainPeaksUncov(i) = sum(NGzPeaksUncov(trainstims))./sum(nbasenoresp(trainstims));
   PGzNovPeaksUncov(i) = sum(NGzPeaksUncov(novstims))./sum(nbasenoresp(novstims));
   %ngztrainuncov2(i)=mean(ngzuncov(trainstims)./nbasenoresp(trainstims));
   %ngznovuncov2(i)=mean(ngzuncov(novstims)./nbasenoresp(novstims));   
   % Omit sites with no baseline bins with no response (2 sites)
   mgztrainuncov(i)=mean(gztrainuncov);
   mgznovuncov(i)=mean(gznovuncov);
   mtgzmot(i)=mean(tgzmot);
   mngzmot(i)=mean(ngzmot);
   tgzmot2=tgzmot;
   ngzmot2=ngzmot;
   tgzmot2(tzeromots)=[];
   ngzmot2(nzeromots)=[];
   mtgzmot2(i)=mean(tgzmot2);% omit motifs with zero FR in both gz and base conditions
   mngzmot2(i)=mean(ngzmot2);
   tbasemot(find(tbasemot==0))=0.1;
   nbasemot(find(nbasemot==0))=0.1;
   tgzmot(find(tgzmot==0))=0.1;
   ngzmot(find(ngzmot==0))=0.1;
   tbasemot2(find(tbasemot2==0))=0.1;
   nbasemot2(find(nbasemot2==0))=0.1;
   tgzmot2(find(tgzmot2==0))=0.1;
   ngzmot2(find(ngzmot2==0))=0.1;
   mtgzpermot(i)=mean(100*(tgzmot-tbasemot)./tbasemot);
   mngzpermot(i)=mean(100*(ngzmot-nbasemot)./nbasemot);
   mtgzpermot2(i)=mean(100*(tgzmot2-tbasemot2)./tbasemot2);%note omitting motifs with no sig response in base or gz condition
   mngzpermot2(i)=mean(100*(ngzmot2-nbasemot2)./nbasemot2);
   ntgzmotuncov(i)=length(find(tgzmot(find(tbasemot<=basethresh))>gzthresh))./length(find(tbasemot<=basethresh));
   mtgzmotuncov(i)=mean(100*(tgzmot(find(tgzmot(find(tbasemot<=basethresh))>gzthresh))-tbasemot(find(tgzmot(find(tbasemot<=basethresh))>gzthresh)))./tbasemot(find(tgzmot(find(tbasemot<=basethresh))>gzthresh)));
   nngzmotuncov(i)=length(find(ngzmot(find(nbasemot<=basethresh))>gzthresh))./length(find(nbasemot<=basethresh));
   mngzmotuncov(i)=mean(100*(ngzmot(find(ngzmot(find(nbasemot<=basethresh))>gzthresh))-nbasemot(find(ngzmot(find(nbasemot<=basethresh))>gzthresh)))./nbasemot(find(ngzmot(find(nbasemot<=basethresh))>gzthresh)));
   MgzMotLifeSparse(i) = LifetimeSparsness(gzmotrate);
   MgzTrainMotLifeSparse(i) = gzTrainMotLifeSparse;
   MgzNovelMotLifeSparse(i) = gzNovelMotLifeSparse;
   gznum=(sum(gzmotrate./(length(gzmotrate)))).^2;
    gzden=sum((gzmotrate.^2)./(length(gzmotrate)));
    mgzlsparsemot(i)=1-(gznum/gzden);
   
   normbasemotrate=100*(basemotrate./max(basemotrate));
   gzmotrateinc=(gzmotrate-basemotrate);
   [r,p ]=corrcoef(basemotrate,gzmotrateinc);
   Rmotgzinc(i)=r(2);
   Pmotgzinc(i)=p(2);
   allbasemot=[allbasemot basemotrate];
   allgzinc=[allgzinc gzmotrate-basemotrate];
   
   mrecper(i)=mean(recper);
   mrecsponper(i)=mean(recsponper);
   mtrainrecper(i)=mean(trainrecper);
   mnovrecper(i)=mean(novrecper);
   mreclsparse(i)=(reclifesparse);
   mzrectrain(i)=mean(zrectrain);
   mzrecnov(i)=mean(zrecnov);
   mpsigrec(i)=psigrecstims;
   mrecbias(i)=recbias;
   mrecrate(i)=mean(recrate);
   mrecsponrate(i)=mean(recsponrate);
   mrecisi(i)=recisi;
   if min(recrate)==0
       recrate(find(recrate==min(recrate)))=0.01;
   end
   mrecdyn(i)=max(recrate)./min(recrate);
   mtrecmot(i)=mean(trecmot);
   mnrecmot(i)=mean(nrecmot);
   recnum=(sum(recmotrate./(length(recmotrate)))).^2;
    recden=sum((recmotrate.^2)./(length(recmotrate)));
    mreclsparsemot(i)=1-(recnum/recden);
    MrecMotLifeSparse(i) = recMotLifeSparse;
   MrecTrainMotLifeSparse(i) = recTrainMotLifeSparse;
   MrecNovelMotLifeSparse(i) = recNovelMotLifeSparse;
   
    if isempty(find(blocktypes==3))==0
        gbindex=[gbindex i];
        mgbper=[mgbper mean(gbper)];
        mgbsponper=[mgbsponper mean(gbsponper)];
        mtraingbper=[mtraingbper mean(traingbper)];
        mnovgbper=[mnovgbper mean(novgbper)];
        mgblsparse=[mgblsparse (gblifesparse)]; 
        mzgbtrain=[mzgbtrain mean(zgbtrain)];
        mzgbnov=[mzgbnov mean(zgbnov)];
        mpsiggb=[mpsiggb psiggbstims];
        mgbbias=[mgbbias gbbias];
        mgbrate=[mgbrate mean(gbrate)];
        [gbp,h]=signrank(baserate,gbrate);
        gbpvals(i)=gbp;
        if gbp<0.05
            if mean(gbrate)>mean(baserate)
                gbinc=[gbinc i];
            end
            if mean(gbrate)<mean(baserate)
                gbdec=[gbdec i];
            end
        end
    end
    
    if isempty(find(blocktypes==4))==0
        sacindex=[sacindex i];
        msacper=[msacper mean(sacper)];
        msacsponper=[msacsponper mean(sacsponper)];
        mtrainsacper=[mtrainsacper mean(trainsacper)];
        mnovsacper=[mnovsacper mean(novsacper)];
        msaclsparse=[msaclsparse (saclifesparse)];
        mzsactrain=[mzsactrain mean(zsactrain)];
        mzsacnov=[mzsacnov mean(zsacnov)];
        mpsigsac=[mpsigsac psigsacstims];
        msacbias=[msacbias sacbias];
        msacrate=[msacrate mean(sacrate)];
        
        [sacp,h]=signrank(baserate,sacrate);
        sacpvals(i)=sacp;
        if sacp<0.05
            if mean(sacrate)>mean(baserate)
                sacinc=[sacinc i];
            end
            if mean(sacrate)<mean(baserate)
                sacdec=[sacdec i];
            end
        end
    end 
    
    if isempty(find(blocktypes==24))==0
        gzsacindex=[gzsacindex i];
        mgzsacper=[mgzsacper mean(gzsacper)];
        mgzsacsponper=[mgzsacsponper mean(gzsacsponper)];
        mtraingzsacper=[mtraingzsacper mean(traingzsacper)];
        mnovgzsacper=[mnovgzsacper mean(novgzsacper)];
        mgzsaclsparse=[mgzsaclsparse (gzsaclifesparse)]; 
        mzgzsactrain=[mzgzsactrain mean(zgzsactrain)];
        mzgzsacnov=[mzgzsacnov mean(zgzsacnov)];
        mpsiggzsac=[mpsiggzsac psiggzsacstims];
        mgzsacbias=[mgzsacbias gzsacbias];
        mgzsacrate=[mgzsacrate mean(gzsacrate)];
        
        [gzsacp1,h]=signrank(baserate,gzsacrate);
        gzsacp1vals(i)=gzsacp1;
        [gzsacp2,h]=signrank(gzrate,gzsacrate);
        gzsacp2vals(i)=gzsacp2;
        if gzsacp2<0.05
            if mean(gzsacrate)>mean(gzrate)
                gzsacinc=[gzsacinc i];
            end
            if mean(gzsacrate)<mean(gzrate)
                gzsacdec=[gzsacdec i];
            end
        end
        
    end
    
    if isempty(find(blocktypes==23))==0
        gzgbindex=[gzgbindex i];
        mgzgbper=[mgzgbper mean(gzgbper)];
        mgzgbsponper=[mgzgbsponper mean(gzgbsponper)];
        mtraingzgbper=[mtraingzgbper mean(traingzgbper)];
        mnovgzgbper=[mnovgzgbper mean(novgzgbper)];
        mgzgblsparse=[mgzgblsparse (gzgblifesparse)]; 
        mzgzgbtrain=[mzgzgbtrain mean(zgzgbtrain)];
        mzgzgbnov=[mzgzgbnov mean(zgzgbnov)];
        mpsiggzgb=[mpsiggzgb psiggzgbstims];
        mgzgbbias=[mgzgbbias gzgbbias];
        mgzgbrate=[mgzgbrate mean(gzgbrate)];
    end
    
    
    
end

ballgzninc=[];
ballgznincindex=[];
for j=1:length(allgzninc)
    ballgzninc=[ballgzninc allgzninc{j}];
    ballgznincindex=[ballgznincindex allgzincindex{j}];
    
    
    
end

nonoresp=find(isnan(ngztrainuncov)|isnan(ngznovuncov));
ngznovuncov(nonoresp)=[];
ngztrainuncov(nonoresp)=[];  
sus=find(strcmp(sitetypes,'su'));
mus=find(strcmp(sitetypes,'mu'));
% uncovsus=sus;
% uncovmus=mus;
% for q=1:length(nonoresp)
%     uncovsus(find(sus==nonoresp(q)))=[];
%     uncovmus(find(mus==nonoresp(q)))=[];
% end

%remove sites with no baseline responses from isi measurements
jo=find(isnan(mbaseisi));
mbaseisi(jo)=[];
mgzisi(jo)=[];
mrecisi(jo)=[];


ventral=find(depths>2580);
dorsal=find(depths<2580);

short=find(strcmp(stimlengths,'s'));
long=find(strcmp(stimlengths,'l'));

sites444=find(strcmp(birdID,'st444'));
sites445=find(strcmp(birdID,'st445'));
sites450=find(strcmp(birdID,'st450'));
sites466=find(strcmp(birdID,'st466'));
sites484=find(strcmp(birdID,'st484'));
sites490=find(strcmp(birdID,'st490'));
sites493=find(strcmp(birdID,'st493'));
sites504=find(strcmpi(birdID,'st504'));
clear sort

%make these for n-way anovas
rectype=[sitetypes'; sitetypes';sitetypes']; % MU or SU
rectype=[rectype; rectype];
tvec=[mzbasetrain';mzgztrain'; mzrectrain'];% train data
nvec=[mzbasenov';mzgznov';mzrecnov']; % nov data
dvec=[tvec;nvec]; % all data
traintype=[ones(length(tvec),1); 1+ones(length(nvec),1)]; % training condition
drugtype=[ones(length(mzbasetrain),1); 1+ones(length(mzbasetrain),1); 2+ones(length(mzbasetrain),1)]; % drug condition


save('BinFRs', 'AllSiteBinFRBase', 'AllSiteBinFRGz')

save('MeansZscores', 'dorsal','ventral','mbaserate', 'mgzrate','mrecrate','mbasetrain','mbasenov','mgztrain','mgznov','mbasesponrate','mgzsponrate','mrecsponrate','MZBaseRate','MZGzRate','MZRecRate','MZBaseTrain2','MZBaseNov2','MZGzTrain2','MZGzNov2','MZBaseSponRate','MZGzSponRate','MZRecSponRate')

save('Selectivity', 'MbaseMotLifeSparse', 'MgzMotLifeSparse', 'MrecMotLifeSparse');

save('SelectivityTraining', 'MbaseTrainMotLifeSparse', 'MbaseNovelMotLifeSparse', 'MgzTrainMotLifeSparse', 'MgzNovelMotLifeSparse','MrecTrainMotLifeSparse', 'MrecNovelMotLifeSparse');
%%
%make for 2way rep measures anova comparing change in FR with Gz for driven
%and spon rates

Y=[mbaserate mgzrate mbasesponrate mgzsponrate]' ;
S=[1:100 1:100 1:100 1:100]' ;
F1=[ones(1,100) ones(1,100)+1 ones(1,100) ones(1,100)+1]' ;
F2=[ones(1,200)+2 ones(1,200)+3]'  ;
VarNames={'Drug'; 'DrivSpon'};
DrivSponStats=rm_anova2(Y,S,F1,F2,VarNames)  ;



figure
subplot(1,4,2)
DataMat=[mbaserate' mgzrate' mbasesponrate' mgzsponrate'];
bar(mean(DataMat),'k')
hold on
errorbar(mean(DataMat),sem(DataMat),'.k')
hold off
set(gca,'xticklabel',{'Base Driven'; 'Gz Driven'; 'Base Spon'; 'Gz Spon'})
ylabel('Firing Rate')


Y=[mbaserate mgzrate mbasesponrate mgzsponrate]' ;
S=[1:100 1:100 1:100 1:100]' ;
F1=[ones(1,100) ones(1,100)+1 ones(1,100) ones(1,100)+1]' ;
F2=[ones(1,200)+2 ones(1,200)+3]'  ;
VarNames={'Drug'; 'DrivSpon'};
DrivSponStats=rm_anova2(Y,S,F1,F2,VarNames)  ;




%%
figure
NormAllSiteBinFRPerIncMat=AllSiteBinFRPerIncMat./BinFRNormMat;
imagesc(BinFRBins,1:100,NormAllSiteBinFRPerIncMat)
colormap(flipud(hot))
axis square
xlabel('% Increase in Firing Rate')
ylabel('Site')







%%

figure
subplot(2,2,1)
hold on
plot(AllSiteBinBaseFRCV(ventral),AllSiteBinGzFRCV(ventral),'ok')
plot(AllSiteBinBaseFRCV(dorsal),AllSiteBinGzFRCV(dorsal),'^k')
xlabel('Inhibition Intact')
ylabel('Inhibition Blocked')
axis square
xlim([0 10])
ylim([ 0 10])
hold on
plot([0 10],[0 10],'k')
hold off

subplot(2,2,3)
hold on
plot(AllSiteBinBaseFRSTD(ventral),AllSiteBinGzFRSTD(ventral),'ok')
plot(AllSiteBinBaseFRSTD(dorsal),AllSiteBinGzFRSTD(dorsal),'^k')
xlabel('STD Inhibition Intact')
ylabel('STD Inhibition Blocked')
axis square
xlim([0 140])
ylim([ 0 140])
hold on
plot([0 140],[0 140],'k')
hold off

figure
hold on
plot([0 150],[0 150],'k')
plot(AllSiteBinBaseFRSTD,AllSiteBinGzFRSTD,'o','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5],'markersize',4)
alpha(0.5)
xlabel('Variance Baseline Firing Rate')
ylabel('Variance Inhibition Blocked Firing Rate')
axis square
xlim([0 150])
ylim([ 0 150])


hold off

%%

VentralCVData=[AllSiteBinBaseTrainFRCV(ventral) AllSiteBinGzTrainFRCV(ventral) AllSiteBinBaseNovFRCV(ventral) AllSiteBinGzNovFRCV(ventral)]' ;
VentralS=[1:40 1:40 1:40 1:40]' ;
DrugCon=[ones(1,40) ones(1,40)+1 ones(1,40) ones(1,40)+1]' ;
TrainCon=[ones(1,80)+2 ones(1,80)+3]'  ;
VarNames={'Drug'; 'Training'};
VentralCVStats=rm_anova2(VentralCVData,VentralS,DrugCon,TrainCon,VarNames)  ;

DorsalCVData=[AllSiteBinBaseTrainFRCV(dorsal) AllSiteBinGzTrainFRCV(dorsal) AllSiteBinBaseNovFRCV(dorsal) AllSiteBinGzNovFRCV(dorsal)]' ;
DorsalS=[1:60 1:60 1:60 1:60]' ;
DrugCon=[ones(1,60) ones(1,60)+1 ones(1,60) ones(1,60)+1]' ;
TrainCon=[ones(1,120)+2 ones(1,120)+3]'  ;
VarNames={'Drug'; 'Training'};
DorsalCVStats=rm_anova2(DorsalCVData,DorsalS,DrugCon,TrainCon,VarNames)  ;

VentralSTDData=[AllSiteBinBaseTrainFRSTD(ventral) AllSiteBinGzTrainFRSTD(ventral) AllSiteBinBaseNovFRSTD(ventral) AllSiteBinGzNovFRSTD(ventral)]' ;
VentralS=[1:40 1:40 1:40 1:40]' ;
DrugCon=[ones(1,40) ones(1,40)+1 ones(1,40) ones(1,40)+1]' ;
TrainCon=[ones(1,80)+2 ones(1,80)+3]'  ;
VarNames={'Drug'; 'Training'};
VentralSTDStats=rm_anova2(VentralSTDData,VentralS,DrugCon,TrainCon,VarNames)  ;

DorsalSTDData=[AllSiteBinBaseTrainFRSTD(dorsal) AllSiteBinGzTrainFRSTD(dorsal) AllSiteBinBaseNovFRSTD(dorsal) AllSiteBinGzNovFRSTD(dorsal)]' ;
DorsalS=[1:60 1:60 1:60 1:60]' ;
DrugCon=[ones(1,60) ones(1,60)+1 ones(1,60) ones(1,60)+1]' ;
TrainCon=[ones(1,120)+2 ones(1,120)+3]'  ;
VarNames={'Drug'; 'Training'};
DorsalSTDStats=rm_anova2(DorsalSTDData,DorsalS,DrugCon,TrainCon,VarNames)  ;


figure
subplot(2,2,1)
hold on;
plot(AllSiteBinBaseTrainFRCV(ventral),AllSiteBinGzTrainFRCV(ventral),'ob')
plot(AllSiteBinBaseNovFRCV(ventral),AllSiteBinGzNovFRCV(ventral),'om')
xlabel('CV Baseline Firing Rate')
ylabel('CV Inhibition Blocked Firing Rate')
axis square
xlim([0 20])
ylim([ 0 20])
plot([0 20],[0 20],'k')
hold off


subplot(2,2,2)
hold on
plot(AllSiteBinBaseTrainFRCV(dorsal),AllSiteBinGzTrainFRCV(dorsal),'ob')
plot(AllSiteBinBaseNovFRCV(dorsal),AllSiteBinGzNovFRCV(dorsal),'om')
xlabel('CV Baseline Firing Rate')
ylabel('CV Inhibition Blocked Firing Rate')
axis square
xlim([0 10])
ylim([ 0 10])
plot([0 10],[0 10],'k')
hold off

subplot(2,2,3)
hold on
plot(AllSiteBinBaseTrainFRSTD(ventral),AllSiteBinGzTrainFRSTD(ventral),'ob')
plot(AllSiteBinBaseNovFRSTD(ventral),AllSiteBinGzNovFRSTD(ventral),'om')
xlabel('STD Inhibition Intact')
ylabel('STD Inhibition Blocked')
axis square
xlim([0 140])
ylim([ 0 140])
plot([0 140],[0 140],'k')
hold off

subplot(2,2,4)
hold on
plot(AllSiteBinBaseTrainFRSTD(dorsal),AllSiteBinGzTrainFRSTD(dorsal),'ob')
plot(AllSiteBinBaseNovFRSTD(dorsal),AllSiteBinGzNovFRSTD(dorsal),'om')
xlabel('STD Inhibition Intact')
ylabel('STD Inhibition Blocked')
axis square
xlim([0 140])
ylim([ 0 140])
plot([0 140],[0 140],'k')
hold off






%%
%Plots Bar Graphs of Zscores of FRs in Baseline Gz and Rec For Ventral and
%Dorsal Sites


figure
subplot(1,4,2)
VentralZRateMat=[MZBaseRate(ventral)' MZGzRate(ventral)' MZRecRate(ventral)'];
bar(mean(VentralZRateMat),'k')
hold on
errorbar(mean(VentralZRateMat),sem(VentralZRateMat),'.k')
hold off
set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'; 'Recovery'})
ylabel('z-score Firing Rate (SDs)')
ylim([-0.6 1])
%axis square

subplot(1,4,4)
DorsalZRateMat=[MZBaseRate(dorsal)' MZGzRate(dorsal)' MZRecRate(dorsal)'];
bar(mean(DorsalZRateMat),'k')
hold on
errorbar(mean(DorsalZRateMat),sem(DorsalZRateMat),'.k')
hold off
set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'; 'Recovery'})
ylabel('z-score Firing Rate (SDs)')
ylim([ -0.6 1])
%axis square

 %Summary plot showing the change in firing rate for each Ventral Site
    %with gabazine
    
    
    VentralRateMatrix = [mbaserate(ventral)' mgzrate(ventral)' mrecrate(ventral)'];
    subplot(1,4,1)
   %hold on
   %bar(mean(VentralRateMatrix),'w','linewidth',2)
   % errorbar(mean(VentralRateMatrix),sem(VentralRateMatrix),'.k','linewidth',2)
    plot([mbaserate(ventral); mgzrate(ventral); mrecrate(ventral)],'-','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5],'markersize',2)
    hold off
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spikes/Sec') 
    VentralRateMat=[mbaserate(ventral)' ones(length(ventral),1) (1:length(ventral))'; mgzrate(ventral)' ones(length(ventral),1)+1 (1:length(ventral))'; mrecrate(ventral)' ones(length(ventral),1)+2 (1:length(ventral))'];
  %axis square
  
  %Summary plot showing the change in firing rate for each Dorsal site
    %with gabazine
    subplot(1,4,3)
    DorsalRateMatrix = [mbaserate(dorsal)' mgzrate(dorsal)' mrecrate(dorsal)'];
    %hold on
    %bar(mean(DorsalRateMatrix),'w','linewidth',2)
    %errorbar(mean(DorsalRateMatrix),sem(DorsalRateMatrix),'.k','linewidth',2)
    plot([mbaserate(dorsal); mgzrate(dorsal); mrecrate(dorsal)],'-','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5],'markersize',2)
    %hold off
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spikes/Sec') 
    
    DorsalRateMat=[mbaserate(dorsal)' ones(length(dorsal),1) (1:length(dorsal))'; mgzrate(dorsal)' ones(length(dorsal),1)+1 (1:length(dorsal))'; mrecrate(dorsal)' ones(length(dorsal),1)+2 (1:length(dorsal))'];
  %axis square

  
  
  
  
  %%
    %Summary plot showing the change in Spontaneous firing rate for each
    %Ventral Site with gabazine
    
    figure
    VentralSponRateMatrix = [mbasesponrate(ventral)' mgzsponrate(ventral)' mrecsponrate(ventral)'];
    subplot(1,4,1)
   % hold on
   %bar(mean(VentralSponRateMatrix),'w','linewidth',2)
   % errorbar(mean(VentralSponRateMatrix),sem(VentralSponRateMatrix),'.k','linewidth',2)
    plot([mbasesponrate(ventral); mgzsponrate(ventral); mrecsponrate(ventral)],'-','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5],'markersize',2)
    %hold off
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spikes/Sec') 
    VentralSponRateMat=[mbasesponrate(ventral)' ones(length(ventral),1) (1:length(ventral))'; mgzsponrate(ventral)' ones(length(ventral),1)+1 (1:length(ventral))'; mrecsponrate(ventral)' ones(length(ventral),1)+2 (1:length(ventral))'];
  
  %Summary plot showing the change in Spontaneous firing rate for each Dorsal site
    %with gabazine
    subplot(1,4,3)
    DorsalSponRateMatrix = [mbasesponrate(dorsal)' mgzsponrate(dorsal)' mrecsponrate(dorsal)'];
    %hold on
    %bar(mean(DorsalSponRateMatrix),'w','linewidth',2)
    %errorbar(mean(DorsalSponRateMatrix),sem(DorsalSponRateMatrix),'.k','linewidth',2)
    plot([mbasesponrate(dorsal); mgzsponrate(dorsal); mrecsponrate(dorsal)],'-','color',[.5 .5 .5],'markerfacecolor',[.5 .5 .5],'markersize',2)
    %hold off
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spikes/Sec') 
    hold off
    DorsalSponRateMat=[mbasesponrate(dorsal)' ones(length(dorsal),1) (1:length(dorsal))'; mgzsponrate(dorsal)' ones(length(dorsal),1)+1 (1:length(dorsal))'; mrecsponrate(dorsal)' ones(length(dorsal),1)+2 (1:length(dorsal))'];    
        
        
%Plots Bar Graphs of Zscores of Spon FRs in Baseline Gz and Rec For Ventral and
%Dorsal Sites


subplot(1,4,2)
VentralZSponRateMat=[MZBaseSponRate(ventral)' MZGzSponRate(ventral)' MZRecSponRate(ventral)'];
bar(mean(VentralZSponRateMat),'k')
hold on
errorbar(mean(VentralZSponRateMat),sem(VentralZSponRateMat),'.k')
hold off
set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'; 'Recovery'})
ylabel('z-score Spontaneous Firing Rate (SDs)')
ylim([-.6 1.0])
%axis square

subplot(1,4,4)
DorsalZSponRateMat=[MZBaseSponRate(dorsal)' MZGzSponRate(dorsal)' MZRecSponRate(dorsal)'];
bar(mean(DorsalZSponRateMat),'k')
hold on
errorbar(mean(DorsalZSponRateMat),sem(DorsalZSponRateMat),'.k')
hold off
set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'; 'Recovery'})
ylabel('z-score SpontaneousFiring Rate (SDs)')
ylim([-.6 1.0])
%ylim([ -0.6 1])
%axis square
%%
VentralBaseDiff = mbasetrain(ventral) - mbasenov(ventral);
VentralGzDiff = mgztrain(ventral) - mgznov(ventral);
VentralBins = linspace(-25,25,30);
VentralBaseHistVals = hist(VentralBaseDiff,VentralBins);
VentralGzHistVals = hist(VentralGzDiff,VentralBins);
VentralDiffMat = [VentralBaseHistVals' -VentralGzHistVals'-VentralBaseHistVals'];
DorsalBaseDiff = mbasetrain(dorsal) - mbasenov(dorsal);
DorsalGzDiff = mgztrain(dorsal) - mgznov(dorsal);
DorsalBins = linspace(-10,10,30);
DorsalBaseHistVals = hist(DorsalBaseDiff,DorsalBins);
DorsalGzHistVals = hist(DorsalGzDiff,DorsalBins);
DorsalDiffMat = [DorsalBaseHistVals' -DorsalGzHistVals'-DorsalBaseHistVals'];
DorsalDiffMat2 = [DorsalGzHistVals' -DorsalBaseHistVals'-DorsalGzHistVals'];

% figure
% subplot(1,2,2)
% bar(VentralBins,VentralDiffMat,'stack');
% xlabel('Learned Firing Rate - Unfamiliar Firing Rate (Spikes/Sec)')
% ylabel('# Sites')
% title('Ventral')
% subplot(1,2,1)
% bar(DorsalBins,DorsalDiffMat2,'stack');
% xlabel('Learned Firing Rate - Unfamiliar Firing Rate (Spikes/Sec)')
% ylabel('# Sites')
% title('Dorsal')
% ylim([-16 16])
% %xlim([-15 15])


% Plots histograms of the difference in FR for learned and unfam songs
% separately for dorsal & ventral and baseline and gabazine
figure
subplot(2,2,1)
hist(DorsalBaseDiff,DorsalBins,'k')
h=findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
set(gca,'tickdir','out')
xlim([-10 10])
ylim([0 16])
hold on
plot([mean(DorsalBaseDiff) mean(DorsalBaseDiff)],[12 14],'k'); 
hold off
xlabel('Increase in Firing Rate to Learned Songs (Spikes/Sec)')
ylabel('# Sites')
title('Dorsal Baseline')
subplot(2,2,3)
hist(DorsalGzDiff,DorsalBins,'k')
h=findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3]);
set(gca,'tickdir','out')
xlim([-10 10])
ylim([0 16])
hold on
plot([mean(DorsalGzDiff) mean(DorsalGzDiff)],[12 14],'k'); 
hold off
xlabel('Increase in Firing Rate to Learned Songs (Spikes/Sec)')
ylabel('# Sites')
title('Dorsal Gabazine')
subplot(2,2,2)
hist(VentralBaseDiff,VentralBins,'k')
h=findobj(gca,'Type','patch');
set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
set(gca,'tickdir','out')
xlim([-25 25])
ylim([0 14])
hold on
plot([mean(VentralBaseDiff) mean(VentralBaseDiff)],[10 12],'k'); 
hold off
xlabel('Increase in Firing Rate to Learned Songs (Spikes/Sec)')
ylabel('# Sites')
title('Ventral Baseline')
subplot(2,2,4)
hist(VentralGzDiff,VentralBins,'k')
h=findobj(gca,'Type','patch');
set(h,'FaceColor',[.8 .3 .3],'EdgeColor',[.8 .3 .3]);
set(gca,'tickdir','out')
xlim([-25 25])
ylim([0 14])
hold on
plot([mean(VentralGzDiff) mean(VentralGzDiff)],[10 12],'k'); 
hold off
xlabel('Increase in Firing Rate to Learned Songs (Spikes/Sec)')
ylabel('# Sites')
title('Ventral Gabazine')

%%

% Scatter Plots of FR to Unfam VS Learn for Base Gz and Diff for Ventral
% and Dorsal Note: Firing Rates have Spontaneous Subtracted
figure
subplot(2,3,1)
plot((mbasenov(ventral)-mbasesponrate(ventral)),(mbasetrain(ventral)-mbasesponrate(ventral)),'ok','markersize',4,'markerfacecolor','k')
xlim([0 160])
ylim([ 0 160])
axis square
hold on
plot([0 160],[0 160],'k')
hold off
xlabel('Firing Rate to Unfamiliar Songs (Spikes/Sec)')
ylabel('Firing Rate to Learned Songs (Spikes/Sec)')

subplot(2,3,2)
plot(mgznov(ventral)-mgzsponrate(ventral),mgztrain(ventral)-mgzsponrate(ventral),'ok','markersize',4,'markerfacecolor','k')
xlim([0 180])
ylim([ 0 180])
axis square
hold on
plot([0 180],[0 180],'k')
hold off
xlabel('Firing Rate to Unfamiliar Songs (Spikes/Sec)')
ylabel('Firing Rate to Learned Songs (Spikes/Sec)')

subplot(2,3,3)
plot((mbasenov(ventral)-mbasesponrate(ventral))-(mbasetrain(ventral)-mbasesponrate(ventral)),(mgznov(ventral)-mgzsponrate(ventral))-(mgztrain(ventral)-mgzsponrate(ventral)),'ok','markersize',4,'markerfacecolor','k')
xlim([-15 20])
ylim([ -15 20])
axis square
hold on
plot([-15 20],[-15 20],'k')
hold off

subplot(2,3,4)
plot((mbasenov(dorsal)-mbasesponrate(dorsal)),(mbasetrain(dorsal)-mbasesponrate(dorsal)),'ok','markersize',4,'markerfacecolor','k')
xlim([0 40])
ylim([ 0 40])
axis square
hold on
plot([0 40],[0 40],'k')
hold off
xlabel('Firing Rate to Unfamiliar Songs (Spikes/Sec)')
ylabel('Firing Rate to Learned Songs (Spikes/Sec)')

subplot(2,3,5)
plot(mgznov(dorsal)-mgzsponrate(dorsal),mgztrain(dorsal)-mgzsponrate(dorsal),'ok','markersize',4,'markerfacecolor','k')
xlim([0 50])
ylim([ 0 50])
axis square
hold on
plot([0 50],[0 50],'k')
hold off
xlabel('Firing Rate to Unfamiliar Songs (Spikes/Sec)')
ylabel('Firing Rate to Learned Songs (Spikes/Sec)')

subplot(2,3,6)
plot((mbasenov(dorsal)-mbasesponrate(dorsal))-(mbasetrain(dorsal)-mbasesponrate(dorsal)),(mgznov(dorsal)-mgzsponrate(dorsal))-(mgztrain(dorsal)-mgzsponrate(dorsal)),'ok','markersize',4,'markerfacecolor','k')
xlim([-10 10])
ylim([ -10 10])
axis square
hold on
plot([-10 10],[-10 10],'k')

hold off

%%
% Plots mean zscore of firing rates with zscore taken across rates from
% both baseline and gabazine conditions
figure
DorsalMZ = [MZBaseTrain2(dorsal)' MZBaseNov2(dorsal)' MZGzTrain2(dorsal)' MZGzNov2(dorsal)'];
VentralMZ = [MZBaseTrain2(ventral)' MZBaseNov2(ventral)' MZGzTrain2(ventral)' MZGzNov2(ventral)'];
subplot(2,1,1)
bar([1 2 5 6],mean(DorsalMZ),'k');
hold on
errorbar([1 2 5 6],mean(DorsalMZ),sem(DorsalMZ),'.k');
xlim([0 7])
set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=60 Dorsal sites')
hold off
axis square
subplot(2,1,2)
bar([1 2 5 6],mean(VentralMZ),'k');
hold on
errorbar([1 2 5 6], mean(VentralMZ),sem(VentralMZ),'.k');
xlim( [0 7])
set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=40 Ventral sites')
hold off
axis square



figure 
DorsalMZ3 = [MZBaseTrain3(dorsal)' MZBaseNov3(dorsal)' MZGzTrain3(dorsal)' MZGzNov3(dorsal)' MZRecTrain3(dorsal)' MZRecNov3(dorsal)'];
VentralMZ3 = [MZBaseTrain3(ventral)' MZBaseNov3(ventral)' MZGzTrain3(ventral)' MZGzNov3(ventral)' MZRecTrain3(ventral)' MZRecNov3(ventral)'];
subplot(2,1,1)
bar([1 2 5 6 9 10],mean(DorsalMZ3),'k');
hold on
errorbar([1 2 5 6 9 10],mean(DorsalMZ3),sem(DorsalMZ3),'.k');
xlim([0 11])
set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar';'Rec-Learned';'Rec-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=60 Dorsal sites')
hold off
axis square
subplot(2,1,2)
bar([1 2 5 6 9 10],mean(VentralMZ3),'k');
hold on
errorbar([1 2 5 6 9 10], mean(VentralMZ3),sem(VentralMZ3),'.k');
xlim( [0 11])
set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'; 'Rec-Learned'; 'Rec-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=40 Ventral sites')
hold off
axis square


%%
% use this stuff to run a repeated measture 2 way anova for dorsal sites
% (rm_anova2)
DorsalZMat=[MZBaseTrain2(dorsal)'; MZBaseNov2(dorsal)'; MZGzTrain2(dorsal)'; MZGzNov2(dorsal)'];
DorsalRMat=[mbasetrain(dorsal)'; mbasenov(dorsal)';mgztrain(dorsal)'; mgznov(dorsal)'];
U=[1:length(dorsal)]';
DorsalSite=[U; U; U; U];
DorsalLearning=[ones(240,1)];
DorsalLearning(61:120)=DorsalLearning(61:120)+1;
DorsalLearning(181:240)=DorsalLearning(181:240)+1;
DorsalDrug=[ones(240,1)];
DorsalDrug(121:240)=DorsalDrug(121:240)+1;
DorsalVarNames={'Learning','Drug Condition'};
Zstats=rm_anova2(DorsalZMat,DorsalSite,DorsalLearning,DorsalDrug,DorsalVarNames)
Rstats=rm_anova2(DorsalRMat,DorsalSite,DorsalLearning,DorsalDrug,DorsalVarNames)

DorsalZMat3=[MZBaseTrain3(dorsal)'; MZBaseNov3(dorsal)'; MZGzTrain3(dorsal)'; MZGzNov3(dorsal)'; MZRecTrain3(dorsal)'; MZRecNov3(dorsal)'];

%%
% use this stuff to run a repeated measture 2 way anova for ventral sites
% (rm_anova2)
VentralZMat=[MZBaseTrain2(ventral)'; MZBaseNov2(ventral)'; MZGzTrain2(ventral)'; MZGzNov2(ventral)'];
VentralRMat=[mbasetrain(ventral)'; mbasenov(ventral)';mgztrain(ventral)'; mgznov(ventral)'];
U=[1:length(ventral)]';
VentralSite=[U; U; U; U];
VentralLearning=[ones(160,1)];
VentralLearning(41:80)=VentralLearning(41:80)+1;
VentralLearning(121:160)=VentralLearning(121:160)+1;
VentralDrug=[ones(160,1)];
VentralDrug(81:160)=VentralDrug(81:160)+1;
VentralVarNames={'Learning','Drug Condition'};
Zstats=rm_anova2(VentralZMat,VentralSite,VentralLearning,VentralDrug,VentralVarNames)
Rstats=rm_anova2(VentralRMat,VentralSite,VentralLearning,VentralDrug,VentralVarNames)
%%

%  Plots mean zscore FR to learned and unfam songs for baseline ang
%  gabazine for dorsal and ventral separately (zcore taken for baseline 
%  and gabazine separately)
figure
subplot(2,1,1)
bar([mean(mzbasetrain(dorsal)) mean(mzbasenov(dorsal)) mean(mzgztrain(dorsal)) mean(mzgznov(dorsal)) ],'k')
    xlim([0 5])
    hold on
    errorbar([mean(mzbasetrain(dorsal)) mean(mzbasenov(dorsal)) mean(mzgztrain(dorsal)) mean(mzgznov(dorsal)) ],[sem(mzbasetrain(dorsal)) sem(mzbasenov(dorsal)) sem(mzgztrain(dorsal)) sem(mzgznov(dorsal)) ],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=60 Dorsal sites')
hold off
subplot(2,1,2)
    bar([mean(mzbasetrain(ventral)) mean(mzbasenov(ventral)) mean(mzgztrain(ventral)) mean(mzgznov(ventral)) ],'k')
    xlim([0 5])
    hold on
    errorbar([mean(mzbasetrain(ventral)) mean(mzbasenov(ventral)) mean(mzgztrain(ventral)) mean(mzgznov(ventral)) ],[sem(mzbasetrain(ventral)) sem(mzbasenov(ventral)) sem(mzgztrain(ventral)) sem(mzgznov(ventral)) ],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('z-score Firing Rate')
    title('n=40 Ventral sites')
    hold off
    
    %%

%  As above for raw firing rates
figure
subplot(2,1,1)
bar([mean(mbasetrain(dorsal)) mean(mbasenov(dorsal)) mean(mgztrain(dorsal)) mean(mgznov(dorsal)) ],'k')
    xlim([0 5])
    hold on
    errorbar([mean(mbasetrain(dorsal)) mean(mbasenov(dorsal)) mean(mgztrain(dorsal)) mean(mgznov(dorsal)) ],[sem(mbasetrain(dorsal)) sem(mbasenov(dorsal)) sem(mgztrain(dorsal)) sem(mgznov(dorsal)) ],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('Firing Rate (Spikes/Sec)')
    title('n=60 Dorsal sites')
hold off
subplot(2,1,2)
    bar([mean(mbasetrain(ventral)) mean(mbasenov(ventral)) mean(mgztrain(ventral)) mean(mgznov(ventral)) ],'k')
    xlim([0 5])
    hold on
    errorbar([mean(mbasetrain(ventral)) mean(mbasenov(ventral)) mean(mgztrain(ventral)) mean(mgznov(ventral)) ],[sem(mbasetrain(ventral)) sem(mbasenov(ventral)) sem(mgztrain(ventral)) sem(mgznov(ventral)) ],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'})
    ylabel('Firing Rate (Spikes/Sec)')
    title('n=40 Ventral sites')
    hold off
%% 
%   Plots mean increase in FR with Gabazine for learned and unfam songs,
%   separately for dorsal and ventral 
VentralDiff=mtraingzdiff(ventral)-mnovgzdiff(ventral);
DorsalDiff=mtraingzdiff(dorsal)-mnovgzdiff(dorsal);
VentralDiffBins=linspace(-10,10,20);
DorsalDiffBins=linspace(-10,10,20);
%     figure
%     subplot(2,2,1)
%     hist(VentralDiff,VentralDiffBins,'k')
%     h=findobj(gca,'Type','patch');
%     set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
%     set(gca,'tickdir','out')
%     xlim([-11 11])
%     ylim([ 0 25])
%     ylabel('Increase in Release from Inhibition to Learned Songs')
%     hold on
%     plot([mean(VentralDiff) mean(VentralDiff)],[10 13],'k'); 
%     hold off
%     
%     subplot(2,2,2)
%     hist(DorsalDiff,DorsalDiffBins,'k')
%     h=findobj(gca,'Type','patch');
%         set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
%     set(gca,'tickdir','out')
%     xlim([-11 11])
%     ylim([0 25])
%     ylabel('Increase in Release from Inhibition to Learned Songs')
%    hold on
%     plot([mean(DorsalDiff) mean(DorsalDiff)],[10 13],'k'); 
%   hold off


figure
subplot(1,2,1)
plot(mnovgzdiff(ventral),mtraingzdiff(ventral),'o')
xlabel('Unfamiliar')
ylabel('Learned')
axis square
xlim([-10 70])
ylim([-10 70])
hold on
plot([-10 70],[-10 70],'k')
hold off

subplot(1,2,2)
plot(mnovgzdiff(dorsal),mtraingzdiff(dorsal),'o')
xlabel('Unfamiliar')
ylabel('Learned')
axis square
xlim([-10 30])
ylim([-10 30])
hold on
plot([-10 30],[-10 30],'k')
hold off




  %%
  
  figure
  VentralTrainGzZ = MTrainGzZDiff(ventral);
  VentralNovGzZ = MNovGzZDiff(ventral);
  DorsalTrainGzZ = MTrainGzZDiff(dorsal);
  DorsalNovGzZ = MNovGzZDiff(dorsal);
  
  
  
  
  
  
  
  
  
  
  
%%

% Same as above but increase in FR plotted as a % of the baseline FR
    figure
    subplot(2,1,1)
    bar([mean(mtraingzper(dorsal)) mean(mnovgzper(dorsal))],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtraingzper(dorsal)) mean(mnovgzper(dorsal)) ],[sem(mtraingzper(dorsal)) sem(mnovgzper(dorsal))],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'})
    ylabel('% Increase in Firing Rate with Gabazine (spikes/sec)')
    title('n=60 Dorsal sites')
    hold off
    subplot(2,1,2)
    bar([mean(mtraingzper(ventral)) mean(mnovgzper(ventral))],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtraingzper(ventral)) mean(mnovgzper(ventral)) ],[sem(mtraingzper(ventral)) sem(mnovgzper(ventral))],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'})
    ylabel('% Increase in Firing Rate with Gabazine (spikes/sec)')
    title('n=40 Ventral sites')

%%

%   Plots ZSCORE of mean increase in FR with Gabazine for learned and unfam songs,
%   separately for dorsal and ventral 
    figure
    subplot(2,1,1)
    bar([mean(MZTrainDiff(dorsal)) mean(MZNovDiff(dorsal))],'k')
    xlim([0 3])
    ylim([-0.3 0.3])
    hold on
    errorbar([mean(MZTrainDiff(dorsal)) mean(MZNovDiff(dorsal)) ],[sem(MZTrainDiff(dorsal)) sem(MZNovDiff(dorsal))],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'})
    ylabel('Z-Score of Increase in Firing Rate with Gabazine (spikes/sec)')
    title('n=60 Dorsal sites')
    hold off
    subplot(2,1,2)
    bar([mean(MZTrainDiff(ventral)) mean(MZNovDiff(ventral))],'k')
    xlim([0 3])
        ylim([-0.3 0.3])
    hold on
    errorbar([mean(MZTrainDiff(ventral)) mean(MZNovDiff(ventral)) ],[sem(MZTrainDiff(ventral)) sem(MZNovDiff(ventral))],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'})
    ylabel('Z-Score of Increase in Firing Rate with Gabazine (spikes/sec)')
    title('n=40 Ventral sites')

%%

%   Plots # of uncovered responses (nonconsecutive peaks) in learned adn
%   unfamiliar responses for dorsal and ventral sites
    figure
    subplot(2,1,1)
    bar([mean(NGzTrainPeaksUncov(dorsal)) mean(NGzNovPeaksUncov(dorsal))],'k')
    xlim([0 3])
    hold on
    errorbar([mean(NGzTrainPeaksUncov(dorsal)) mean(NGzNovPeaksUncov(dorsal)) ],[sem(NGzTrainPeaksUncov(dorsal)) sem(NGzNovPeaksUncov(dorsal))],'.k')
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('# of Uncovered Responses')
    title('n=60 Dorsal sites')
    hold off
    subplot(2,1,2)
    bar([mean(NGzTrainPeaksUncov(ventral)) mean(NGzNovPeaksUncov(ventral))],'k')
    xlim([0 3])
    hold on
    errorbar([mean(NGzTrainPeaksUncov(ventral)) mean(NGzNovPeaksUncov(ventral)) ],[sem(NGzTrainPeaksUncov(ventral)) sem(NGzNovPeaksUncov(ventral))],'.k')
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('# of Uncovered Responses')
    title('n=40 Ventral sites')

    %%

%   Plots Frequency of Peaks with uncovered responses (nonconsecutive peaks) in learned adn
%   unfamiliar responses for dorsal and ventral sites
    figure
    subplot(2,2,1)
    DorsalFGzBins=linspace(-1.2,1.2,30);
    DorsalFGzDiff=FGzTrainPeaks(dorsal)-FGzNovPeaks(dorsal);
    hist(DorsalFGzDiff,DorsalFGzBins)
    h=findobj(gca,'Type','patch');
    set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
    set(gca,'tickdir','out')
    xlabel('Frequency of Uncovered Responses (Peaks/Sec)')
    xlim([-1.25 1.25])
    hold on
    plot([mean(DorsalFGzDiff) mean(DorsalFGzDiff)],[10 13],'k'); 
    hold off
    
    subplot(2,2,2)
    VentralFGzBins=linspace(-.5,.5,30);
    VentralFGzDiff=FGzTrainPeaks(ventral)-FGzNovPeaks(ventral);
    hist(VentralFGzDiff,VentralFGzBins)
    h=findobj(gca,'Type','patch');
    set(h,'FaceColor',[.5 .5 .5],'EdgeColor',[.5 .5 .5]);
    set(gca,'tickdir','out')
    xlim([-.6 .6])
    xlabel('Frequency of Uncovered Responses (Peaks/Sec)')
    hold on
    plot([mean(VentralFGzDiff) mean(VentralFGzDiff)],[10 13],'k'); 
    hold off
%%
    
    figure
    subplot(1,2,2)
    DorsalFreqPeaksMat=[FGzTrainPeaks(dorsal)' FGzNovPeaks(dorsal)'];
    VentralFreqPeaksMat=[FGzTrainPeaks(ventral)' FGzNovPeaks(ventral)'];    
    bar(mean(DorsalFreqPeaksMat),'k');
    hold on;
    errorbar(mean(DorsalFreqPeaksMat),sem(DorsalFreqPeaksMat),'.k')
    hold off
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Frequency of Uncovered Responses')
    
    subplot(1,2,1)
    bar(mean(VentralFreqPeaksMat),'k');
    hold on;
    errorbar(mean(VentralFreqPeaksMat),sem(VentralFreqPeaksMat),'.k')
    hold off
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Frequency of Uncovered Responses')
    
    
    
    
    figure
    subplot(1,4,2)
    DorsalNumPeaksMat=[NGzTrainPeaksUncov(dorsal)' NGzNovPeaksUncov(dorsal)'];
    VentralNumPeaksMat=[NGzTrainPeaksUncov(ventral)' NGzNovPeaksUncov(ventral)'];    
    bar(mean(DorsalNumPeaksMat),'k');
    hold on;
    errorbar(mean(DorsalNumPeaksMat),sem(DorsalNumPeaksMat),'.k')
    hold off
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Number of Responses')
    
    subplot(1,4,1)
    bar(mean(VentralNumPeaksMat),'k');
    hold on;
    errorbar(mean(VentralNumPeaksMat),sem(VentralNumPeaksMat),'.k')
    hold off
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Number of Responses')
%%





% Scatterplot of the FR of each motif in the baseline condition vs. the
% increae in FR with gabazine
figure
plot(allbasemot,allgzinc,'ok')
 xlabel('Baseline Motif Firing Rate spikes/sec')
 ylabel('Increase in Firing Rate with Gabazine spikes/sec')
 




%%
% Shows the range of the % the increase in FR with gabazine
% when respsones are binned
    [x,y]=sort((AllSiteGzIncMax-AllSiteGzIncMin)./2);
     figure
    subplot(2,2,1);
    MinMaxGzIncMat=[AllSiteGzIncMin(y)' AllSiteGzIncMax(y)'];
     YLabMat=[(1:100)' (1:100)'];
     plot(MinMaxGzIncMat',YLabMat','-ok','markerfacecolor','k','markersize',2)
     
     ylabel('Sites')
     xlabel('Increase in Firing Rate with Gabazine (Spikes/Sec)')
%     
        subplot(2,2,2)
        [x,y]=sort((AllSiteGzPerIncMax-AllSiteGzPerIncMin)./2);
    MinMaxGzPerIncMat=[AllSiteGzPerIncMin(y)' AllSiteGzPerIncMax(y)'];
     YLabMat=[(1:100)' (1:100)'];
     plot(MinMaxGzPerIncMat',YLabMat','-ok','markerfacecolor','k','markersize',2)
     
     ylabel('Sites')
     xlabel('%Increase in Firing Rate with Gabazine (Spikes/Sec)')
    
    
    
    
    
 %%   
    
    %  As above for the raw increase in FR with gabazine
    [x,y]=sort(medgzinc);
    figure
    errorbar(medgzinc(y),iqrgzinc(y)./2,'ok','markerfacecolor','k');
    %set(gcf,'markerfacecolor','k')
    xlabel('Sites')
    ylabel('Change in Firing Rate with Gabazine (median+iqr)')
    
    
  
    
    %%
  
  % Plots the lifetime sparsenesss for baseline, gabazine and recovery
  % condition
    figure
    plot([mbaselsparse; mgzlsparse; mreclsparse],'-o')
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Lifetime Sparseness') 
    title('Single/Multi Units n=75')
    hold on
    errorbar([mean(mbaselsparse) mean(mgzlsparse) mean(mreclsparse)],[sem(mbaselsparse) sem(mgzlsparse) sem(mreclsparse)],'-k','linewidth',6)
  hold off  
    
  %%
  
  % Plots the dynamic range for baseline, gabazine and recovery condition,
  % not as useful as lifetime sparseness
    figure
    bar([mean(mbasedyn) mean(mgzdyn) mean(mrecdyn)],'k')
    hold on
    errorbar([mean(mbasedyn) mean(mgzdyn) mean(mrecdyn)],[sem(mbasedyn) sem(mgzdyn) sem(mrecdyn)],'.k')
  hold off
    xlim([0 4])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Dynamic Range in Firing Rate') 
    title('Multi/Single Units n=75')
    
    
    
    
    
  %%  
    figure
    % Summary plot showing hte change in FR of each site with GABA added
    
    %gbsus=intersect(sus,gbindex);
    %gbmus=inersect(mus,gbindex);
    plot([mbaserate(gbindex); mgzrate(gbindex); mgbrate; mrecrate(gbindex)],'-o')
    xlim([0.5 4.5])
    set(gca,'xtick',[1 2 3 4])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Gaba'; 'Recovery'})
    ylabel('Mean Firing Rate') 
    title('Single/Multi Units n=16')
    hold on
    errorbar([mean(mbaserate(gbindex)) mean(mgzrate(gbindex)) mean(mgbrate) mean(mrecrate(gbindex))],[sem(mbaserate(gbindex)) sem(mgzrate(gbindex)) sem(gbrate) sem(mrecrate(gbindex))],'-k','linewidth',6)
  hold off
    GbRateMat=[mbaserate(gbindex)' ones(length(gbindex),1) (1:length(gbindex))'; mgbrate' ones(length(gbindex),1)+1 (1:length(gbindex))'; mrecrate(gbindex)' ones(length(gbindex),1)+2 (1:length(gbindex))'];
       
       
    %%
  
  % Summary plot showing the change in FR of each site with Saclofen added
    figure
    plot([mbaserate(sacindex); mgzrate(sacindex); msacrate; mrecrate(sacindex)],'-o')
    xlim([0.5 4.5])
    set(gca,'xtick',[1 2 3 4])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Saclofen'; 'Recovery'})
    ylabel('Mean Firing Rate') 
    title('Single/Multi Units n=6')
    hold on
    errorbar([mean(mbaserate(sacindex)) mean(mgzrate(sacindex)) mean(msacrate) mean(mrecrate(sacindex))],[sem(mbaserate(sacindex)) sem(mgzrate(sacindex)) sem(msacrate) sem(mrecrate(sacindex))],'-k','linewidth',6)
  hold off
   
  SacRateMat=[mbaserate(sacindex)' ones(length(sacindex),1) (1:length(sacindex))'; msacrate' ones(length(sacindex),1)+1 (1:length(sacindex))'; mrecrate(sacindex)' ones(length(sacindex),1)+2 (1:length(sacindex))'];
    
    
    
  
  %%
Y=[MbaseTrainMotLifeSparse MbaseNovelMotLifeSparse MgzTrainMotLifeSparse MgzNovelMotLifeSparse]' ;
S=[1:100 1:100 1:100 1:100]' ;
F1=[ones(1,100) ones(1,100)+1 ones(1,100) ones(1,100)+1]' ;
F2=[ones(1,200)+2 ones(1,200)+3]'  ;
VarNames={'Training'; 'Drug'};
SelectivityTrainingStats=rm_anova2(Y,S,F1,F2,VarNames)  ;
  
Y=[MbaseTrainMotLifeSparse(dorsal) MbaseNovelMotLifeSparse(dorsal) MgzTrainMotLifeSparse(dorsal) MgzNovelMotLifeSparse(dorsal)]' ;
S=[1:60 1:60 1:60 1:60]' ;
F1=[ones(1,60) ones(1,60)+1 ones(1,60) ones(1,60)+1]' ;
F2=[ones(1,120)+2 ones(1,120)+3]'  ;
VarNames={'Training'; 'Drug'};
DorsalSelectivityTrainingStats=rm_anova2(Y,S,F1,F2,VarNames)  ;

Y=[MbaseTrainMotLifeSparse(ventral) MbaseNovelMotLifeSparse(ventral) MgzTrainMotLifeSparse(ventral) MgzNovelMotLifeSparse(ventral)]' ;
S=[1:40 1:40 1:40 1:40]' ;
F1=[ones(1,40) ones(1,40)+1 ones(1,40) ones(1,40)+1]' ;
F2=[ones(1,80)+2 ones(1,80)+3]'  ;
VarNames={'Training'; 'Drug'};
VentralSelectivityTrainingStats=rm_anova2(Y,S,F1,F2,VarNames)  ;


    figure
    subplot(1,4,1)
    bar([mean(MbaseMotLifeSparse) mean(MgzMotLifeSparse) mean(MrecMotLifeSparse)],'k')
   % xlim([0.5 3.5])
   % set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'; 'Recovery'})
    ylabel('Motif Lifetime Sparseness') 
   % title('Single/Multi Units n=75')
    hold on
    errorbar([mean(MbaseMotLifeSparse) mean(MgzMotLifeSparse) mean(MrecMotLifeSparse)],[sem(MbaseMotLifeSparse) sem(MgzMotLifeSparse) sem(MrecMotLifeSparse)],'.k')
  hold off
  
  subplot(1,4,2)
    bar([mean(MbaseMotLifeSparse) mean(MgzMotLifeSparse) ],'k')
   % xlim([0.5 3.5])
   % set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'})
    ylabel('Motif Lifetime Sparseness') 
   % title('Single/Multi Units n=75')
    hold on
    errorbar([mean(MbaseMotLifeSparse) mean(MgzMotLifeSparse) ],[sem(MbaseMotLifeSparse) sem(MgzMotLifeSparse) ],'.k')
  hold off
  
  subplot(1,4,3)
  bar([mean(MbaseMotLifeSparse(ventral)) mean(MgzMotLifeSparse(ventral)) ],'k')
  set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'})
  ylabel('Motif Lifetime Sparseness') 
  hold on
  errorbar([mean(MbaseMotLifeSparse(ventral)) mean(MgzMotLifeSparse(ventral)) ],[sem(MbaseMotLifeSparse(ventral)) sem(MgzMotLifeSparse(ventral)) ],'.k')
  ylim([0 0.5])
  hold off
  
  subplot(1,4,4)
  bar([mean(MbaseMotLifeSparse(dorsal)) mean(MgzMotLifeSparse(dorsal)) ],'k')
  set(gca,'xticklabel',{'Baseline'; 'Inhibition Blocked'})
  ylabel('Motif Lifetime Sparseness') 
  hold on
  errorbar([mean(MbaseMotLifeSparse(dorsal)) mean(MgzMotLifeSparse(dorsal)) ],[sem(MbaseMotLifeSparse(dorsal)) sem(MgzMotLifeSparse(dorsal)) ],'.k')
  hold off
  
  %%
  figure
  learnMotLSparse=[MbaseTrainMotLifeSparse' MbaseNovelMotLifeSparse' MgzTrainMotLifeSparse' MgzNovelMotLifeSparse' MrecTrainMotLifeSparse' MrecNovelMotLifeSparse'];
    bar(mean(learnMotLSparse),'k')
    xlim([0 7])
    set(gca,'xtick',[1 2 3 4 5 6])
    set(gca,'xticklabel',{'Train Baseline'; 'Novel Baseline'; 'Train Gabazine'; 'Novel Gabazine'; 'Train Recovery'; 'Novel Recovery'})
    ylabel('Motif Lifetime Sparseness') 
    title('Single/Multi Units n=75')
    hold on
    errorbar(mean(learnMotLSparse), sem(learnMotLSparse),'.k')
  hold off
    %%
    figure
    plot([mbaselsparsemot; mgzlsparsemot; mreclsparsemot],'-o')
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Motif Lifetime Sparseness') 
    title('Single/Multi Units n=75')
    hold on
    errorbar([mean(mbaselsparsemot) mean(mgzlsparsemot) mean(mreclsparsemot)],[sem(mbaselsparsemot) sem(mgzlsparsemot) sem(mreclsparsemot)],'-k','linewidth',6)
  hold off
    
    
    %%
  
  % Summary plot showing the change in FR for each site when Gabazine and
  % Saclofen are added together
    figure
    plot([mbaserate(gzsacindex); mgzrate(gzsacindex); mgzsacrate; mrecrate(gzsacindex)],'-o')
    xlim([0.5 4.5])
    set(gca,'xtick',[1 2 3 4])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Gabazine+Saclofen'; 'Recovery'})
    ylabel('Mean Firing Rate') 
    title('Single/Multi Units n=13')
    hold on
    errorbar([mean(mbaserate(gzsacindex)) mean(mgzrate(gzsacindex)) mean(mgzsacrate) mean(mrecrate(gzsacindex))],[sem(mbaserate(gzsacindex)) sem(mgzrate(gzsacindex)) sem(mgzsacrate) sem(mrecrate(gzsacindex))],'-k','linewidth',6)
  hold off
  GzSacRateMat=[mbaserate(gzsacindex)' ones(length(gzsacindex),1) (1:length(gzsacindex))'; mgzrate(gzsacindex)' ones(length(gzsacindex),1)+1 (1:length(gzsacindex))'; mgzsacrate' ones(length(gzsacindex),1)+2 (1:length(gzsacindex))'; mrecrate(gzsacindex)' ones(length(gzsacindex),1)+3 (1:length(gzsacindex))'];
    
    
%%

%   Plots the Spontaneous Firing Rate of each site
    figure
    plot([mbasesponrate; mgzsponrate; mrecsponrate],'-o')
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spontaneous Firing Rate') 
    title('Single/Multi Units n=75')
    
    hold on
        errorbar([mean(mbasesponrate) mean(mgzsponrate) mean(mrecsponrate)],[sem(mbasesponrate) sem(mgzsponrate) sem(mrecsponrate)],'-k','linewidth',6)
        hold off
   
        
  %%
       
    
     figure
    errorbar([mean(mbasesponrate) mean(mgzsponrate) mean(mrecsponrate)],[sem(mbasesponrate) sem(mgzsponrate) sem(mrecsponrate)],'-xk','linewidth',4)
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Spontaneous Firing Rate (Spikes/Sec)') 
    title('Single/Multi Units n=75')  
    
       
    %%
    altmbasesponrate=mbasesponrate;
    altmbasesponrate(find(altmbasesponrate==0))=1;
    normmbasesponrate=100*(altmbasesponrate./altmbasesponrate);
    normmgzsponrate=100*(mgzsponrate./altmbasesponrate);
    normmrecsponrate=100*(mrecsponrate./altmbasesponrate);
    
    
     figure
    plot([normmbasesponrate; normmgzsponrate; normmrecsponrate],'-o')
    xlim([0.5 3.5])
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('% of Baseline Spontaneous Firing Rate') 
    title('Single/Multi Units n=75')
 
    
    
    
    
    
    
 
 %%   
    figure
    bar([mean(mtraingzdiff) mean(mnovgzdiff)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtraingzdiff) mean(mnovgzdiff)],[sem(mtraingzdiff) sem(mnovgzdiff)],'.k')
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Change in Firing Rate with Gabazine (spikes/sec)')
    title('Single/Multi Units n=100')
    
    
 
    
    
    
    
   
 %%   
    figure
    bar([mean(ngztrainuncov) mean(ngznovuncov)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(ngztrainuncov) mean(ngznovuncov)],[sem(ngztrainuncov) sem(ngznovuncov)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Proportion of Uncovered Response') 
    title('Multi/Single Units n=75')
%% 
  figure
  rid=find(isnan(mgztrainuncov));
  mgztrainuncov(rid)=[];
  mgznovuncov(rid)=[];
  
    bar([mean(mgztrainuncov) mean(mgznovuncov)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mgztrainuncov) mean(mgznovuncov)],[sem(mgztrainuncov) sem(mgznovuncov)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Mean Firing Rate of Uncovered Response (spikes/sec)') 
    title('Multi/Single Units n=69')



    
    
    
    
    
    
    



%%

vsus=intersect(ventral,sus);
figure
    bar([mean(mzbasetrain(vsus)) mean(mzbasenov(vsus)) mean(mzgztrain(vsus)) mean(mzgznov(vsus)) mean(mzrectrain(vsus)) mean(mzrecnov(vsus))],'k')
    xlim([0 7])
    hold on
    errorbar([mean(mzbasetrain(vsus)) mean(mzbasenov(vsus)) mean(mzgztrain(vsus)) mean(mzgznov(vsus)) mean(mzrectrain(vsus)) mean(mzrecnov(vsus))],[sem(mzbasetrain(vsus)) sem(mzbasenov(vsus)) sem(mzgztrain(vsus)) sem(mzgznov(vsus)) sem(mzrectrain(vsus)) sem(mzrecnov(vsus))],'.k')
    set(gca,'xticklabel',{'Base-Learned'; 'Base-Unfamiliar'; 'Gz-Learned'; 'Gz-Unfamiliar'; 'Rec-Learned'; 'Rec-Unfamilliar'})    
    ylabel('z-score Firing Rate Baseline Condition')
    title('n=5 SU ventral')
%%


figure
    bar([mean(mtgzmot-mtbasemot) mean(mngzmot-mnbasemot)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtgzmot-mtbasemot) mean(mngzmot-mnbasemot)],[sem(mtgzmot-mtbasemot) sem(mngzmot-mnbasemot)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Increase in Motif Firing Rate with Gabazine spikes/sec') 
    title('Multi/Single Units n=75')
    
    %%
    
    figure
    bar([mean(100*(mtgzmot-mtbasemot)./mtbasemot) mean(100*(mngzmot-mnbasemot)./mnbasemot)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(100*(mtgzmot-mtbasemot)./mtbasemot) mean(100*(mngzmot-mnbasemot)./mnbasemot)],[sem(100*(mtgzmot-mtbasemot)./mtbasemot) sem(100*(mngzmot-mnbasemot)./mnbasemot)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('%Increase in Motif Firing Rate with Gabazine spikes/sec') 
    title('Multi/Single Units n=75')
    
    
    
    
    
    
    
    %%
rems=isnan(mngzmot2)|isnan(mtgzmot2)|isnan(mnbasemot2)|isnan(mtbasemot2);
mtgzmot2(rems)=[];
mtbasemot2(rems)=[];
mngzmot2(rems)=[];
mnbasemot2(rems)=[];

figure
    bar([mean(mtgzmot2-mtbasemot2) mean(mngzmot2-mnbasemot2)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtgzmot2-mtbasemot2) mean(mngzmot2-mnbasemot2)],[sem(mtgzmot2-mtbasemot2) sem(mngzmot2-mnbasemot2)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('Increase in Firing Rate with Gabazine spikes/sec') 
    title('Multi/Single Units n=73')

%%
rems2= union(find(isnan(mngzpermot2)),find(isnan(mtgzpermot2)));
mtgzpermot2(rems2)=[];
mngzpermot2(rems2)=[];



figure
    bar([mean(mtgzpermot2) mean(mngzpermot2)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtgzpermot2) mean(mngzpermot2)],[sem(mtgzpermot2) sem(mngzpermot2)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('% Increase in Firing Rate with Gabazine spikes/sec') 
    title('Multi/Single Units n=73')
    
    %%
    
   figure
    bar([mean(mtgzpermot) mean(mngzpermot)],'k')
    xlim([0 3])
    hold on
    errorbar([mean(mtgzpermot) mean(mngzpermot)],[sem(mtgzpermot) sem(mngzpermot)],'.k')
    set(gca,'xtick',[1 2 ])
    set(gca,'xticklabel',{'Learned'; 'Unfamiliar'})
    ylabel('% Increase in Motif Firing Rate with Gabazine spikes/sec') 
    title('Multi/Single Units n=75') 

    
    %%   
    figure
    bar([mean(mbaseisi) mean(mgzisi) mean(mrecisi)],'k')
    xlim([0 4])
    hold on
    errorbar([mean(mbaseisi) mean(mgzisi) mean(mrecisi)],[sem(mbaseisi) sem(mgzisi) sem(mrecisi)],'.k')
    set(gca,'xtick',[1 2 3])
    set(gca,'xticklabel',{'Baseline'; 'Gabazine'; 'Recovery'})
    ylabel('Mean ISI') 
    title('Multi/Single Units n=74')
%%
posbias=find(mbasebias>0);
negbias=find(mbasebias<0);

figure
    bar([mean(mzbasetrain(posbias)) mean(mzbasenov(posbias)) mean(mzgztrain(posbias)) mean(mzgznov(posbias)) mean(mzrectrain(posbias)) mean(mzrecnov(posbias)) mean(mzbasetrain(negbias)) mean(mzbasenov(negbias)) mean(mzgztrain(negbias)) mean(mzgznov(negbias)) mean(mzrectrain(negbias)) mean(mzrecnov(negbias))],'k')
   % xlim([0 5])
    hold on
    errorbar([mean(mzbasetrain(posbias)) mean(mzbasenov(posbias)) mean(mzgztrain(posbias)) mean(mzgznov(posbias)) mean(mzrectrain(posbias)) mean(mzrecnov(posbias)) mean(mzbasetrain(negbias)) mean(mzbasenov(negbias)) mean(mzgztrain(negbias)) mean(mzgznov(negbias)) mean(mzrectrain(negbias)) mean(mzrecnov(negbias))],[sem(mzbasetrain(posbias)) sem(mzbasenov(posbias)) sem(mzgztrain(posbias)) sem(mzgznov(posbias)) sem(mzrectrain(posbias)) sem(mzrecnov(posbias)) sem(mzbasetrain(negbias)) sem(mzbasenov(negbias)) sem(mzgztrain(negbias)) sem(mzgznov(negbias)) sem(mzrectrain(negbias)) sem(mzrecnov(negbias))],'.k')
    set(gca,'xticklabel',{'Base-Learned +bias'; 'Base-Unfamiliar+bias'; 'GZ-Learned+bias'; 'GZ-Unfamiliar+bias'; 'Rec-Learned+bias'; 'Rec-Unfamiliar+bias'; 'Base-Learned -bias'; 'Base-Unfamiliar-bias'; 'GZ-Learned-bias'; 'GZ-Unfamiliar-bias'; 'Rec-Learned-bias'; 'Rec-Unfamiliar-bias'})
    ylabel('z-score Firing Rate ')
   % title('sites preferring unfamiliar songs n=36')


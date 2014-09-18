
% [wavname, wav, FS] = batchGetStimsionto

allnames_ionto = [];
allwavs_ionto = [];

% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1230/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1460/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1590/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1810/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2570/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2680/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z3200/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2000/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2450/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1550/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,0);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];
% what are the "start and stop" and blocksize is 10... but there are 30
% trials?
% this is probably because there are 10 trials for each of 3 periods of
% pharma?
% [wavname, wav, FS] = batchGetStimsionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1670/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1750/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1860/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2150/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2330/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2770/ss001
[wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1000/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1055/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1151/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1260/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1330/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1440/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1480/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1570/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1700/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z2650/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% % cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z780/ss001
% % [wavname, wav, FS] = batchGetStimsionto2
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z920/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,5,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss003
% [wavname, wav, FS] = batchGetStimsionto(-2,12,5,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss004
% [wavname, wav, FS] = batchGetStimsionto(-2,12,5,[1 2 5],'su',0,1)
% 
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1110/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1150/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1260/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1310/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1440/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1520/ss003
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1640/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1760/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2190/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2250/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1370/ss001
% [wavname, wav, FS] = batchGetStimsionto(-2,12,10,[1 2 5],'mu',0,1)


cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z830/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2040/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2100/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2250/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2350/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2480/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 3 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2600/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2800/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];



cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z780/ss002/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 24 5 3 23 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2150/ss002/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 24 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2210/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2250/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2311/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 24 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2430/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2800/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5 4 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2880/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];



cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss002/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2150/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2280/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2430/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2650/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2810/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2910/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2990/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z3032/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];



cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z970/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1170/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1461/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1570/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1850/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2101/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2220/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 5],'su',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2620/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2811/ss001/
[wavname, wav, FS] = batchGetStimsionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
allnames_ionto = [allnames_ionto, wavname];
allwavs_ionto = [allwavs_ionto, wav];











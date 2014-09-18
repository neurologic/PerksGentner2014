%%%%%%%%%
% re-saving 'batchbatchPSTHionto' as 'SpikingResponse_InhibitionEffect' 
% and cleaning up the plotting and analyses parts
%%%%%%%%%%%%%

load([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
%%
load([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_ResultBaseRef.mat'],'ResultBaseRef')
load([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_ResultGzRef.mat'],'ResultGzRef')


%% alert tone
[alert,fsalert] = wavread('/Users/kperks/Dropbox/presentations/FigsForRecruitment/media2.wav');
alert = alert(end-(2*fsalert):end,1);
%%
bps = 10000;
sigon = 20000;
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1230/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1460/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1590/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1810/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2570/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2680/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z3200/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2000/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2450/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
base_data_All = [];
gz_data_All = [];

isite = 1
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1550/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,0);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];
% [base_data, gz_data] = batchPSTHionto(bps,start,stop,blocksize,blocktype,sitetype,smooth,doprintout)

%2
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1670/ss001
[base_data2, gz_data2] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%3
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1750/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];
%this one seems to become more bursty with removal of inhibition...
%and reduction in spike rates? at least majorly shrunken windows...
%complete transformation to bursty responses...
%can do a check for this? ... that only analyze sites with the same isi
%distribution between cases?

%it still has gain >1 , but it also seems to asymptote? not good linear?

%4
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1860/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%5
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2150/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%6
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2330/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%7
cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2770/ss001
[base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%sound(alert,fsalert)
save([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')

% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1000/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1055/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1151/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1260/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1330/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1440/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1480/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1570/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1700/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z2650/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% % cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z780/ss001
% % batchrasterionto2
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z920/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,5,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss003
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,5,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss004
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,5,[1 2 5],'su',0,1)
%
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1110/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1150/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1260/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1310/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1440/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1520/ss003
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1640/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1760/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2190/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2250/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1370/ss001
% [base_data, gz_data] = batchPSTHionto(bps,-2,12,10,[1 2 5],'mu',0,1)

%8
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z830/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%9
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2040/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%9
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2100/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%10
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2250/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%11
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2350/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%12
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2480/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 3 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%13
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2600/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%14
cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2800/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


save([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
%sound(alert,fsalert)

%15
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z780/ss002/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 24 5 3 23 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%16
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2150/ss002/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 24 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%17
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2210/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%18
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2250/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%19
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2311/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 24 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%20
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2430/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%21
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2800/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5 4 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%22
cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2880/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

save([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
%sound(alert,fsalert)



%23
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%24
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss002/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%25
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2150/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%26
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2280/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%27
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2430/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%28
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2650/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%29
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2810/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%30
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2910/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%31
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2990/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%32
cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z3032/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

save([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
%sound(alert,fsalert)



%33
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z970/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%34
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1170/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%35
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1461/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%36
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1570/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


%%%%%%%%%%pick up here... wed april 2nd when get to campus
load([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
isite = 37
%37
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1850/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%38
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2101/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%39
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2220/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 5],'su',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%40
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2620/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%41
cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2811/ss001/
[base_data, gz_data] = batchPSTHionto(bps,-2,7,10,[1 2 24 5 3 5],'mu',0,1);
Result(isite) = SpkResponseGain(base_data,gz_data,sigon,bps);
isite = isite +1
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


save([r.Dir.Expt 'Analysis/Meta_Responses_DB/SpikesGainAnal_Results.mat'],'Result')
%sound(alert,fsalert)
%sound(alert,fsalert)

%% plot the cdf of all responses in each condition and baseline rate in each condition as hist

allintactrate_resp =  [];
allgzrate_resp =  [];
allintactrate_pre = [];
allgzrate_pre = [];

for isite = 2:size(Result,2) %%%site 2 seems to be a repeat of 1
    allintactrate_resp =  [allintactrate_resp,Result(isite).baseRate_resp];
    allgzrate_resp =  [allgzrate_resp,Result(isite).gzRate_resp];
    
    allintactrate_pre = [allintactrate_pre,Result(isite).baseRate_pre];
    allgzrate_pre = [allgzrate_pre,Result(isite).gzRate_pre];
end

figure;
hold on
[x,p] = empcdf(allintactrate_resp);
stairs(x,p,'color','k')
[x,p] = empcdf(allgzrate_resp);
stairs(x,p,'color','r')
xlabel('response spike rate')
legend({'intact','gz'})
ylabel('cdf')
medianRate_stim_intact = median(allintactrate_resp)
CI = getCDFconf(allintactrate_resp,95)
medianRate_stim_gz = median(allgzrate_resp)
CI = getCDFconf(allgzrate_resp,95)
nresponses = size(allintactrate_resp)
nresponses = size(allgzrate_resp)
kstest(allgzrate_resp)

figure;
hold on
scatter(repmat(1,size(allintactrate_pre,1),size(allintactrate_pre,2)),allintactrate_pre)
scatter(repmat(2,size(allgzrate_pre,1),size(allgzrate_pre,2)),allgzrate_pre,'MarkerEdgeColor','r')

D = [allintactrate_pre',allgzrate_pre'];
hfig = figure;
boxplot(D,'factorgap',10,'color','kr')
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/BoxPlotPreStimRates_intactVgz.jpg')
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/BoxPlotPreStimRates_intactVgz.fig')
medianRate_pre_intact = median(allintactrate_pre)
CI = getCDFconf(allintactrate_pre,95)
medianRate_pre_gz = median(allgzrate_pre)
CI = getCDFconf(allgzrate_pre,95)
nstimuli = size(allintactrate_pre)
nstimuli = size(allgzrate_pre)
%%
% doplot =1;
% if doplot ==1
%     prefig = figure;hold on
%     title('pre-stim')
%     line([0,50],[0,50],'color','k')
%     stimfig = figure;hold on
%     title('stim')
%     line([0,250],[0,250],'color','k')
%     gainfig = figure;hold on
%     title('correlation between gain and magnitude of response at baseline')
% end
% 
% mn_resp_gain_dec = [];
% mn_resp_scale_dec = [];
% mn_baseresp_dec = [];
% mn_resp_gain_inc = [];
% mn_resp_base = [];
% mn_resp_gain = [];
% fit_resp_gain = [];
% mn_pre_gain = [];


p_corr = nan(size(Result,2),1);
R2 = nan(size(Result,2),1);
IO_scale = nan(size(Result,2),1);
IO_offset = nan(size(Result,2),1);
fits = nan(size(Result,2),2); %(column 1 is multiplicative; column 2 is subtract)
fits_CI(1).CI = nan(size(Result,2),2); %with field CI: array(1) is mult; array(2) is subt
%column 1 is lower CI bound, column 2 is upper CI bound
fits_CI(2).CI = nan(size(Result,2),2);
signif_mult = nan(size(Result,2),1); %1 if significant slope of fit
signif_sub = nan(size(Result,2),1); %0 if significant offset of fit
resid_slopeConstant = [];
resid_offsetConstant = [];

slopefig = figure;hold on;title('fitlines with significant slopes')
set(gca,'XLim',[0,2],'YLim',[0,2]);
line([0,2],[0,2],'color','k')
offsetfig = figure;hold on;title('fitlines with significant offsets')
set(gca,'XLim',[0,2],'YLim',[0,2]);
line([0,2],[0,2],'color','k')

scat = figure; hold on
set(gca,'XLim',[0,2],'YLim',[0,2]);
line([0,2],[0,2],'color','k')

scales = [];
diffs = [];

fig_sortscat = figure; hold on
fig_IOfits = figure; hold on
%%%%%%before run this make sure that Result is 40 long (take away the first
%%%%%%one because first and second sites are repeats
for isite = 2:size(Result,2) %%%site 2 seems to be a repeat of 1
    allintactrate =  Result(isite).baseRate_resp;
    allgzrate =  Result(isite).gzRate_resp;
    
    respinds = find(allintactrate>0);
    allintactrate = allintactrate(respinds);
    allgzrate = allgzrate(respinds);
    
    if size(allintactrate,2) > 10
    scale =  max(allgzrate);
    allgzrate = allgzrate ./ scale;
    allintactrate = allintactrate ./ scale;
    figure(scat);scatter(allintactrate,allgzrate)
    %     if doplot ==1
    %         figure(stimfig)
    %         scatter(allbaserate,allgzrate)
    %     end
    
    
        %     options = fitoptions('poly1');
        
        [r,p] = corrcoef(allintactrate',allgzrate');
        p_corr(isite) = p(~triu(p));
        R2(isite) = (r(~triu(r)))^2;
        cf = fit(allintactrate',allgzrate','poly1');
        fitline = cf.p1*(allintactrate) + cf.p2;
        resid_fit{isite} = allgzrate - fitline;
        CI = confint(cf,0.95);
        %slope CI is column 1
        %offset CI is column 2
        fits (isite,1) = cf.p1;%(column 1 is multiplicative; column 2 is subtract)
        fits (isite,2) = cf.p2;
        fits_CI(1).CI(isite,:) = CI(:,1)';
        fits_CI(2).CI(isite,:) = CI(:,2)';
        scales = [scales,(allgzrate./allintactrate)];
        diffs = [diffs,(allgzrate - allintactrate)];
        

        if isempty(crossing(CI(:,1),[],1)) %if the CI are both on either side of 1
            signif_mult(isite,1) = 1; %1 if significant slope of fit
                        figure(slopefig)
                        xaxis = [0,1];
%                         line(allbaserate,fitline)
                        if cf.p1 > 1
                        line(xaxis, cf.p1*(xaxis) + cf.p2,'color','r');
                        end
                        if cf.p1 < 1
                        line(xaxis, cf.p1*(xaxis) + cf.p2,'color','b');
                        end
        elseif ~isempty(crossing(CI(:,1),[],1)) %slope not different from 1
            signif_mult(isite,1) = 0;
        end
        if isempty(crossing(CI(:,2),[],0)) %if the CI are both on either side of 1
            signif_sub(isite,1) = 1; %1 if significant slope of fit
                        figure(offsetfig)
                        xaxis = [0,1];
%                         line(allbaserate,fitline)
                        if cf.p1 > 1
                        line(xaxis, cf.p1*(xaxis) + cf.p2,'color','r');
                        end
                        if cf.p1 < 1
                        line(xaxis, cf.p1*(xaxis) + cf.p2,'color','b');
                        end
        elseif ~isempty(crossing(CI(:,2),[],0)) %slope not different from 1
            signif_sub(isite,1) = 0;
        end
        
        
        %get a fit for the sorted responses based on gz rates
        %this is like a pseudo I/O curve
%         for each site, plot the scatter
        [b,ix] = sort(allintactrate);
        xaxis = [1:size(ix,2)];
        figure(fig_sortscat)
        subplot(7,7,isite); hold on
        scatter(xaxis,allintactrate(ix))
        cf = fit(xaxis',allintactrate(ix)','poly1');
        baseIO = cf.p1*(xaxis) + cf.p2;
        scatter(xaxis,allgzrate(ix))
        line(xaxis,baseIO,'color','k','LineWidth',3)
        cf = fit(xaxis',allgzrate(ix)','poly1');
        gzIO = cf.p1*(xaxis) + cf.p2;
        line(xaxis,gzIO,'color','r','LineWidth',3)
        
        figure(fig_IOfits)
        scatter(baseIO,gzIO)
        cf = fit(baseIO',gzIO','poly1');
        IO_scale(isite) = cf.p1;
        IO_offset(isite) = cf.p2;
%         line(xaxis,baseIO,'color','k')
%         line(xaxis,gzIO,'color','r')
        
        options = fitoptions('poly1');
        %force fit to slope == 1; calculate residuals
        options.Upper = [1,Inf];
        options.Lower = [1,-Inf];
        cf = fit(allintactrate',allgzrate','poly1',options);
        fitline = cf.p1*(allintactrate) + cf.p2;
        resid_slopeConstant{isite} = allgzrate - fitline;
        
        
        %force fit to offset == 0; calculate residuals
        options.Upper = [Inf,0];
        options.Lower = [-Inf,0];
        cf = fit(allintactrate',allgzrate','poly1',options);
        fitline = cf.p1*(allintactrate) + cf.p2;
        resid_offsetConstant{isite} = allgzrate - fitline;
    end
end

hfig = figure;hold on
set(gca,'XLim',[0.5,3.5],'XTick',[1,2,3],'XTickLabel',{'offset = 0','slope = 1','fit'})
ylabel('mean abs(residuals) of linear fit')
for isite = 1:size(Result,2)
    if ~isempty(resid_offsetConstant{isite}) 
    %remind myself why this would not have been calculated for some sites...
        scatter([1,2],[mean(abs(resid_offsetConstant{isite})),mean(abs(resid_slopeConstant{isite}))],50,'k','fill')
        if kstest2(abs(resid_offsetConstant{isite}),abs(resid_slopeConstant{isite})) == 1
            line([1,2],[mean(abs(resid_offsetConstant{isite})),mean(abs(resid_slopeConstant{isite}))],'color','k')
        end
        scatter(3,mean(abs(resid_fit{isite})),50,'k','fill')
    end
end
% saveas(hfig, '/Volumes/kperks/Documents/EJN/Figures/EPSforFig4/GainModels_Residuals.fig')
% saveas(hfig, '/Volumes/kperks/Documents/EJN/Figures/EPSforFig4/GainModels_Residuals.jpg')
% saveas(hfig, '/Volumes/kperks/Documents/EJN/Figures/EPSforFig4/GainModels_Residuals.eps')

figure;
% c = [min(scales):1:max(scales)];
c = [0:0.2:10];
n = hist(scales,c);
stairs(c,n./size(scales,2),'color','g');
line([1,1],[0,max(n./size(scales,2))],'color','k')
xlabel('gabazine / intact')
ylabel('p')
figure;
c = [min(diffs):0.1:max(diffs)];
n = hist(diffs,c);
stairs(c,n./size(diffs,2),'color','b');
xlabel('gabazine - intact')
ylabel('p')


slfits = fits(find(signif_mult==1),1);
offits = fits(find(signif_sub==1),2);
figure; hold on
c = [0:.2:12]
% c = [min([min(slfits),min(IO_scale)]):0.5:max([max(slfits),max(IO_scale)])];
n = hist(slfits,c);
stairs(c,n./size(slfits,1),'color','k');
c = [0:0.5:12];
n = hist(IO_scale(~isnan(IO_scale)),c);
stairs(c,n./size(find(~isnan(IO_scale)),1),'color','r');
legend({'data','IO fitlines'})
xlabel('slopes from fits')
ylabel('p(slope)')
set(gca,'XTick',[0,1,5,10])

figure; hold on
c = [-3:0.1:3];
% c = [min([min(offits),min(IO_offset)]):0.1:max([max(offits),max(IO_offset)])];
n = hist(offits,c);
stairs(c,n./size(offits,1),'color','k');
n = hist(IO_offset(~isnan(IO_offset)),c);
stairs(c,n./size(find(~isnan(IO_offset)),1),'color','r');
legend({'data','IO fitlines'})
xlabel('offsets from fits')
ylabel('p(offset)')

figure;hold on
scatter(repmat(1,size(offits,1),size(offits,2)),offits)
scatter(repmat(2,size(slfits,1),size(slfits,2)),slfits)
title('significant values from fits')
set(gca,'XLim',[0.5,2.5],'XTick',[1,2],'XTickLabel',...
    {['offsets n=' num2str(size(offits,1))],['slopes n=' num2str(size(slfits,1))] })

fitlogical = ~isnan(p_corr);
fitinds = find(fitlogical);
num_sites_fit = size(fitinds,1)
signif_corr_inds = find(p_corr<=0.05);
%number of significant fits (determined by corr_coeff p val)
signif_fit_inds = intersect(fitinds,signif_corr_inds);
num_signif_fits = size(signif_corr_inds,1)
%-->of those fits that are significant,
%1)how many sites had significant offset
sigoffinds = find(signif_sub(signif_corr_inds));
num_signif_offset = size(sigoffinds,1)
%2)how many sites had significant slope
sigscaleinds = find(signif_mult(signif_corr_inds));
num_signif_scale = size(sigscaleinds,1)
%3)how many sites had both significant slope and offset
signifbothinds = intersect(sigoffinds,sigscaleinds);
num_signif_both = size(signifbothinds,1)


%     options = fitoptions('poly1');
%     options.Upper = [Inf,Inf];
%     options.Lower = [-Inf,-Inf];
%     if ~isempty(Result(isite).baseRate_resp)
%         cf = fit(Result(isite).baseRate_resp',allgzrate','poly1',options);
%         fit_resp_gain (isite) = cf.p1;
%     elseif isempty(Result(isite).baseRate_resp)
%         fit_resp_gain (isite) = NaN;
%     end
%
%     mn_resp_gain (isite) = mean(Result(isite).gzRate_resp ./ Result(isite).baseRate_resp);
%     incinds = find((Result(isite).gzRate_resp-Result(isite).baseRate_resp) > 0);
%     baserespinds = find(Result(isite).baseRate_resp ~= 0);
%     gainincinds = intersect(incinds,baserespinds);
%     mn_resp_gain_inc  = [mn_resp_gain_inc,(Result(isite).gzRate_resp (incinds) ...
%         ./ Result(isite).baseRate_resp(incinds))];
%     decinds = find((Result(isite).gzRate_resp-Result(isite).baseRate_resp) < 0);
%     mn_resp_gain_dec  = [mn_resp_gain_dec,(Result(isite).gzRate_resp (decinds) ...
%         ./ Result(isite).baseRate_resp(decinds))];
%     mn_resp_scale_dec = [mn_resp_scale_dec,(Result(isite).gzRate_resp (decinds) ...
%         - Result(isite).baseRate_resp(decinds))];
%     mn_baseresp_dec = [mn_baseresp_dec,Result(isite).baseRate_resp(decinds)];
%
%     if doplot ==1
%         figure(gainfig)
%     end
%     mn_resp_base (isite) = mean(Result(isite).baseRate_resp);
%     if doplot ==1
%         scatter(mn_resp_base (isite) ,mn_resp_gain (isite))
%     end
%
%     if doplot ==1
%         figure(prefig)
%         scatter(Result(isite).baseRate_pre ,Result(isite).gzRate_pre)
%         cf = fit(Result(isite).baseRate_pre',Result(isite).gzRate_pre','poly1');
%     end
%     mn_pre_gain (isite) = nanmean((Result(isite).gzRate_pre) ./ (Result(isite).baseRate_pre));
%     allprebase = [allprebase, (Result(isite).baseRate_pre)];
%     allpregz = [allpregz, (Result(isite).gzRate_pre)];
% end

%%
%relationship between gain of response and magnitude
respinds = find(allintactrate ~= 0 );
figure
scatter(allintactrate(respinds),(allgzrate(respinds) ./ allintactrate(respinds)))
ylabel('gain')
xlabel('baseline response rate')
allgain = (allgzrate(respinds) ./ allintactrate(respinds));
maxgain = max(allgain)
meangain = mean(allgain)
medgain = median(allgain)

figure
scatter(allintactrate(respinds),allgzrate(respinds))

%scatter and distribution of subtraction
figure;
scatter(allintactrate, (allgzrate - allintactrate))

figure;
edges = [-200:20:200]+10;
% hist((allgzrate - allbaserate),edges)
n = histc((allgzrate - allintactrate),edges)
bar(edges,n./size(allgzrate,2))

%distribution of gains
edges = [0:0.2:10];
n = hist(allgain,edges);
figure;
bar(edges,(n/size(allgain,2)))


%distribution of gains during pre-stim
respinds = find(allprebase ~= 0 );
allpregain = allpregz(respinds) ./ allprebase(respinds);
edges = [0:0.2:10];
n = hist(allpregain,edges);
figure;
bar(edges,(n/size(allpregain,2)))


%% distribution of response strengths
figure 
hold all
siteind = 1;
for isite = 1:size(ResultBaseRef,2)
    if size(ResultBaseRef(isite).baseRate_resp,2) >30
        isnorm(siteind) = kstest(ResultBaseRef(isite).baseRate_resp);
        siteind = siteind + 1;
    [x,p] = empcdf(ResultBaseRef(isite).baseRate_resp);
    stairs(x,p)
    
    end
end

figure 
hold all
for isite = 1:size(ResultBaseRef,2)
    if size(ResultBaseRef(isite).gzRate_resp,2) >30
    [x,p] = empcdf(ResultBaseRef(isite).gzRate_resp);
    stairs(x,p)
    
    end
end

%ranked responses in each condition compared
figure;
hold all
for isite = 1:size(ResultBaseRef,2)
    if size(ResultBaseRef(isite).baseRate_resp,2) >30
       [B,IXbase] = sort(ResultBaseRef(isite).baseRate_resp);
        [B,IXgz] = sort(ResultBaseRef(isite).gzRate_resp);
        scatter(IXbase,IXgz)
      [RHO(isite),PVAL(isite)] = corr(IXbase',IXgz','type','Spearman')  ;
%     [x,p] = empcdf(ResultBaseRef(isite).baseRate_resp);
%     stairs(x,p)
    
    end
end


figure 
hold all
siteind = 1;
for isite = 1:size(ResultGzRef,2)
    if size(ResultGzRef(isite).baseRate_resp,2) >30
        isnorm(siteind) = kstest(ResultGzRef(isite).baseRate_resp);
        siteind = siteind + 1;
    [x,p] = empcdf(ResultGzRef(isite).baseRate_resp);
    stairs(x,p)
    
    end
end

figure 
hold all
for isite = 1:size(ResultGzRef,2)
    if size(ResultGzRef(isite).gzRate_resp,2) >50
    [x,p] = empcdf(ResultGzRef(isite).gzRate_resp);
    stairs(x,p)
    
    end
end

%ranked responses in each condition compared
figure;
hold all
for isite = 1:size(ResultGzRef,2)
    if size(ResultGzRef(isite).baseRate_resp,2) >30
       [B,IXbase] = sort(ResultGzRef(isite).baseRate_resp);
        [B,IXgz] = sort(ResultGzRef(isite).gzRate_resp);
        scatter(IXbase,IXgz)
      [RHO(isite),PVAL(isite)] = corr(IXbase',IXgz','type','Spearman')  ;
%     [x,p] = empcdf(ResultBaseRef(isite).baseRate_resp);
%     stairs(x,p)
    
    end
end

%% distribution of response duratinos

%% distribution of spike times within each response
respind = 1;
for isite = 1:size(Result,2)
    for iresp = 1:size(Result(isite).basespkT,2)
        MnB(respind) = mean(Result(isite).basespkT{iresp});
        MnG(respind) = mean(Result(isite).gzspkT{iresp});
        MnShift(respind) = mean(Result(isite).gzspkT{iresp}) - mean(Result(isite).basespkT{iresp});
        MedShift(respind) = median(Result(isite).gzspkT{iresp}) - median(Result(isite).basespkT{iresp});
        SkewB(respind) = skewness(Result(isite).basespkT{iresp});
        SkewG(respind) = skewness(Result(isite).gzspkT{iresp});
        Skewshift(respind) = SkewG(iresp) - SkewB(iresp);
        respind = respind +1;
    end
    
end

%% isi distribution during intact vs during gz need this during spontaneous

%% fit between gain or diff and magnitude
magR = [];
gainR = [];
diffR = [];
allResp_gz = [];
allResp_base = [];

for isite = 1:size(Result,2)
    magR{isite} = Result(isite).gzRate_resp;
    diffR{isite} = Result(isite).gzRate_resp - Result(isite).baseRate_resp;
    gainR{isite} = Result(isite).gzRate_resp ./ Result(isite).baseRate_resp;
       
    allResp_gz = [allResp_gz, Result(isite).gzRate_resp];
    allResp_base = [allResp_base, Result(isite).baseRate_resp];
end

minmax = [min([allResp_gz,allResp_base]),max([allResp_gz,allResp_base])];
edges = [minmax(1):10:minmax(2)];
figure;
hold on
n = histc(allResp_gz,edges);
stairs(edges,n/size(allResp_gz,2),'color','r','LineWidth',4)
n = histc(allResp_base,edges);
stairs(edges,n/size(allResp_base,2),'color','b','LineWidth',4)

figure;
hold on
[x,p] = empcdf(allResp_gz)
stairs(x,p,'color','r','LineWidth',4)
[x,p] = empcdf(allResp_base)
stairs(x,p,'color','b','LineWidth',4)
xlabel('spike rate','FontSize',14)
ylabel('cummulative fraction','FontSize',14)
legend({'gabazine','baseline'},'FontSize',14)

%%
numresps = 10;
slfig = [];
ll = [];
M_G = [];
M_D = [];
ciM_G = [];
ciM_D = [];
O_G = [];
O_D = [];

prefig = figure;
hold on
stimfig = figure('name',['minreps for fits = ' num2str(numresps)])
hold on
for isite = 1:size(Result,2)
    figure(stimfig)
    slfig(isite) = subplot(6,7,isite);
    hold on
    scatter(magR{isite},gainR{isite},30,'r')
    scatter(magR{isite},diffR{isite},30,'g')
    axis tight
    
    naninf = intersect(find(~isnan(gainR{isite})), find(~isinf(gainR{isite})));
    gainRfilt = gainR{isite}(naninf);
    magRfilt = magR{isite}(naninf);
    if size(naninf,2)<numresps
        ciM_G(isite,:) = [nan,nan];
        M_G(isite) = nan;
        O_G(isite) = nan;
    end
    if size(naninf,2)>=numresps
        fitresult = fit(magRfilt',gainRfilt','poly1');
        yfit = (fitresult.p1 .* magRfilt) + fitresult.p2;
        ci = confint(fitresult);
        ciM_G(isite,:) = ci(:,1);
        M_G(isite) = fitresult.p1;
        ll(1) = line(magRfilt,yfit,'color','r');
        
        sigresult = ((ciM_G(isite,1) > 0 ) && (ciM_G(isite,2) > 0)) || ((ciM_G(isite,1) < 0 ) && (ciM_G(isite,2) < 0)) ;
        %sigresult is 0 if the ci values span across zero
        if  ~sigresult %if the slope is not different from zero, what is that constant value?
            O_G(isite) = fitresult.p2 ;
        end
    end
    
    
    naninf = intersect(find(~isnan(diffR{isite})), find(~isinf(diffR{isite})));
    diffRfilt = diffR{isite}(naninf);
    magRfilt = magR{isite}(naninf);
    if size(naninf,2)<numresps
        ciM_D(isite,:) = [nan,nan];
        M_D(isite) = nan;
        O_D(isite) = nan;
    end
    if size(naninf,2)>=numresps
        fitresult = fit(magRfilt',diffRfilt','poly1');
        yfit = (fitresult.p1 .* magRfilt) + fitresult.p2;
        ci = confint(fitresult);
        ciM_D(isite,:) = ci(:,1);
        M_D(isite) = fitresult.p1;
        ll(2)= line(magRfilt,yfit,'color','g');
        
        sigresult = ((ciM_D(isite,1) > 0 ) && (ciM_D(isite,2) > 0)) || ((ciM_D(isite,1) < 0 ) && (ciM_D(isite,2) < 0)) ;
       
        %sigresult is 0 if the ci values span across zero
        if  ~sigresult %if the slope is not different from zero, what is that constant value?
            O_D(isite) = fitresult.p2 ;
        end
    end
    
    title(['G: [' num2str(ciM_G(isite,:)) ']  D: [' ...
        num2str(ciM_D(isite,:)) ']'], 'FontSize',14)
    
   
end

%CDF of slopes for fit between Gain or Diff and magnitude of response
figure('color','w')
hold on
[x,p]=empcdf(M_G(~isnan(M_G)))
stairs(x,p,'color','m','LineWidth',4)
[x,p]=empcdf(M_D(~isnan(M_D)))
stairs(x,p,'color','g','LineWidth',4)
title(['minreps for fits = ' num2str(numresps)])
ylabel('cummulative fraction','FontSize',14)
xlabel('slope of fit','FontSize',14)
legend({'gain vs magnitude','difference vs magnitude'},'FontSize',14)

%histogram of slopes for fit between Gain or Diff and magnitude of response
figure;
hold on
minmax = [floor(min([min(M_G),min(M_D)])), ceil(max([max(M_G),max(M_D)]))];
edges = [minmax(1):0.2:minmax(2)];
n = histc(M_G(~isnan(M_G)),edges);
stairs(edges,n/size(M_G(~isnan(M_G)),2),'color','m','LineWidth',4)
n = histc(M_D(~isnan(M_D)),edges);
stairs(edges,n/size(M_D(~isnan(M_D)),2),'color','g','LineWidth',4)
title(['minreps for fits = ' num2str(numresps)])
legend({'how much gain varies with mag', 'how much diff varies with mag'},'FontSize',14)
xlabel('fit slope','FontSize',14)
ylabel('normalized frequency','FontSize',14)


%cell array of results
% number of sites with enough reps,
% number of fits with non-zero slope, 
% mean offset when non-zero slope,
% 

fitGcell{1,1} = 'Gain';
fitGcell{2,1} = ['n site > ' num2str(numresps) 'reps'];
nsites = size(M_G,2) - size(find(isnan(M_G)),2);
fitGcell{2,2} = nsites;
fitGcell{3,1} = ['n slope signif 0'];
numslopes = nsites - size(find(O_G ==0),2);
fitGcell{3,2} = numslopes;
fitGcell{4,1} = 'mean slope';
fitGcell{4,2} = nanmean(M_G);
fitGcell{5,1} = 'median slope';
fitGcell{5,2} = median(M_G(~isnan(M_G)));
fitGcell{6,1} = 'mean gain when m=0';
meanslope = mean(O_G(intersect(find(~isnan(O_G)),find(O_G ~=0))));
fitGcell{6,2} = meanslope;
fitGcell{7,1} = 'median gain when m=0';
medslope = median(O_G(intersect(find(~isnan(O_G)),find(O_G ~=0))));
fitGcell{7,2} = medslope;
fitGcell

fitDcell{1,1} = 'Difference';
fitDcell{2,1} = ['n site > ' num2str(numresps) 'reps'];
nsites = size(M_D,2) - size(find(isnan(M_D)),2);
fitDcell{2,2} = nsites;
fitDcell{3,1} = ['n slope signif 0'];
numslopes = nsites - size(find(O_D ==0),2);
fitDcell{3,2} = numslopes;
fitDcell{4,1} = 'mean slope';
fitDcell{4,2} = nanmean(M_D);
fitDcell{5,1} = 'median slope';
fitDcell{5,2} = median(M_D(~isnan(M_D)));
fitDcell{6,1} = 'mean diff when m=0';
meanslope = mean(O_D(intersect(find(~isnan(O_D)),find(O_D ~=0))));
fitDcell{6,2} = meanslope;
fitDcell{7,1} = 'median diff when m=0';
medslope = median(O_D(intersect(find(~isnan(O_D)),find(O_D ~=0))));
fitDcell{7,2} = medslope;
fitDcell


figpre = figure;
hold on
diffP = [];
gainP = [];
for isite = 1:size(Result,2)
 figure(figpre)
  slpre(isite) = subplot(6,7,isite);
%  hold on
 
 maxrate = max([max(Result(isite).gzRate_pre),max(Result(isite).baseRate_pre)]);
 normgz = Result(isite).gzRate_pre ./ maxrate;
 normbase = Result(isite).baseRate_pre ./ maxrate;
%  scatter(normbase,normgz,30)
 
 scatter(Result(isite).gzRate_pre,Result(isite).baseRate_pre,30,'k')
 diffP(isite) = mean(Result(isite).gzRate_pre - Result(isite).baseRate_pre);
gainP(isite) = mean(Result(isite).gzRate_pre ./ Result(isite).baseRate_pre);
magP(isite) = mean(Result(isite).gzRate_pre);
end

a = find(~isinf(gainP));
b = find(~isnan(gainP));
figure;hist(gainP(intersect(a,b)))
figure;scatter(magP(intersect(a,b)),gainP(intersect(a,b)))
cf = fit(magP(intersect(a,b))',gainP(intersect(a,b))','poly1')
ciM = confint(cf)
MG = cf.p1
OG = cf.p2
title(['ci on slope = ' num2str(ciM(:,1)')] )
text(25,6,['gain = ' num2str(cf.p2) '      ci on gain estimate = ' num2str(ciM(:,2)')])

a = find(~isinf(diffP));
b = find(~isnan(diffP));
figure;hist(diffP(intersect(a,b)))
figure;scatter(magP(intersect(a,b)),diffP(intersect(a,b)))
cf = fit(magP(intersect(a,b))',diffP(intersect(a,b))','poly1')
ciD = confint(cf)
MD = cf.p1
OD = cf.p2
title(['ci on slope = ' num2str(ciD(:,1)')] )
text(5,30,['diff = ' num2str(cf.p2) '      ci on diff estimate = ' num2str(ciD(:,2)')])


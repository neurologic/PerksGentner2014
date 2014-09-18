%%%%%%%%%
% re-saving 'batchbatchFANOionto' as 'SpikingReproducibility' 
% and cleaning up the plotting and analyses parts
%%%%%%%%%%%%%

% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1230/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1460/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1590/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z1810/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2570/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z2680/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-300_L1000/Z3200/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2000/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B444/Left_R-500_L1000/Z2450/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
base_data_All = [];
gz_data_All = [];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1550/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,0);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];
% [base_data, gz_data] = batchFANOionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1670/ss001
[base_data2, gz_data2] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1750/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z1860/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2150/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2330/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B445/Left_R-200_L1000/Z2770/ss001
[base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1000/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1055/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1151/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1260/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1330/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1440/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1480/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1570/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z1700/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z2650/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% % cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z780/ss001
% % batchrasterionto2
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z920/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,5,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss003
% [base_data, gz_data] = batchFANOionto(-2,12,5,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B493/Right_R-350_L1100/Z950/ss004
% [base_data, gz_data] = batchFANOionto(-2,12,5,[1 2 5],'su',0,1)
% 
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1110/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1150/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1260/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1310/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1440/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1520/ss003
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1640/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1760/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2190/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z2250/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'su',0,1)
% cd /Users/kperks/GitHub/iontoncm/B484/Left_R-300_L1000/Z1370/ss001
% [base_data, gz_data] = batchFANOionto(-2,12,10,[1 2 5],'mu',0,1)


cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z830/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2040/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2100/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2250/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2350/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2480/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 3 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2600/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B450/Right_R-300_L1000/Z2800/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z780/ss002/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 24 5 3 23 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2150/ss002/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 24 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2210/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2250/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2311/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 24 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-200_L1000/Z2430/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2800/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5 4 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B466/Left_R-350_L1000/Z2880/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z660/ss002/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2150/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2280/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2430/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2650/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2810/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2910/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z2990/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B504/LEFT_R-300_L1000/Z3032/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];


cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z970/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1170/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1461/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1570/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z1850/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2101/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2220/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 5],'su',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2620/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

cd /Users/kperks/GitHub/iontoncm/B490/Right_R-300_L1000/Z2811/ss001/
[base_data, gz_data] = batchFANOionto(-2,7,10,[1 2 24 5 3 5],'mu',0,1);
base_data_All = [base_data_All base_data];
gz_data_All = [gz_data_All gz_data];

%%
meanspks = [];
tmpspks = [];
trial_len = [];
spkvec = [];
for icond = 1:max(size(base_data_All))
  tmpspks = zeros(size(base_data_All(icond).spikes,1),size(base_data_All(icond).spikes,2));
  for itrial = 1:size(base_data_All(icond).spikes,1)
    tmpspks(itrial,find(base_data_All(icond).spikes(itrial,:))) = 1;  
  end
  trial_len(icond) = size(base_data_All(icond).spikes,2);
  meanspks{icond} = mean(tmpspks);
end

for icond = 1:max(size(base_data_All))
    spkvec(icond,:) = meanspks{icond}(1,1:min(trial_len));
end

%% set up params

times  = 25:50:size(tmp_pre,2)-24; % from 1950 ms before target onset until 2000 ms after.
fanoParams.alignTime = size(tmp_pre,2)/2; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
fanoParams.matchReps = 10;
fanoParams.binSpacing = 0.25;
fanoParams.boxWidth = boxWidth;

%% split and combine gz and base data
% to get pre and post stim array that has 
% base as first half and gz for second half of split for Fano anal
start_stim = 2000;
boxWidth = 50;
data_len = (boxWidth * floor((start_stim - 1)/50));
pre_data = [];
post_data = [];
tmp_post = [];
tmp_pre = [];
for iar = 1:size(gz_data_All,2)
    this_ar_gz = gz_data_All(iar).spikes;
    this_ar_bs = base_data_All(iar).spikes;
    
    for itrial = 1:size(this_ar_gz,1)
        prestart = (start_stim - 1) - data_len;
    tmp_pre(itrial,:) = [this_ar_bs(itrial,prestart:(start_stim - 1)) ...
        this_ar_gz(itrial,prestart:(start_stim - 1)) ];
    
        postend = (start_stim + 50) + data_len;
    tmp_post(itrial,:) = [this_ar_bs(itrial,(start_stim + 50):4000) ...
        this_ar_gz(itrial,(start_stim + 50):4000) ];
    
    end 
    
    pre_data(iar).spikes = tmp_pre;
    post_data(iar).spikes = tmp_post;
    
end
%%

%%%%%%%%%%%%%%%
%for gaba condition...
times  = 25:50:size(tmp_pre,2)-24; % from 1950 ms before target onset until 2000 ms after.
fanoParams.alignTime = (size(tmp_pre,2)/2)-25; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
fanoParams.matchReps = 10;
fanoParams.binSpacing = 0.25;
fanoParams.boxWidth = boxWidth;
Result_pre = VarVsMean(pre_data, times, fanoParams); 
plotFano(Result_pre);

hfig = figure;hold on
subplot(4,1,1)
hold on
data_split = size(Result_pre.meanRateAll,1)/2;

plot(Result_pre.meanRateSelect(1:data_split),'color','k','LineWidth',3)
plot(Result_pre.meanRateSelect(data_split + 1 :end),'color','r','LineWidth',3)
axis tight
ylabel('mean count Raw','FontSize',16)

subplot(4,1,2);hold on
plot(Result_pre.FanoFactor(1:data_split),'color','k','LineWidth',3)
plot(Result_pre.FanoFactor(data_split + 1 :end),'color','r','LineWidth',3)
axis tight
ylabel('FanoFactor matched','FontSize',16)

subplot(4,1,3);hold on
range = max(Result_pre.FanoFactor)-min(Result_pre.FanoFactor);
bins = (range / 10);
edges = [min(Result_pre.FanoFactor):bins:max(Result_pre.FanoFactor)];
bins_b = histc(Result_pre.FanoFactor(1:data_split),edges)./data_split;
bins_i = histc(Result_pre.FanoFactor(data_split + 1 :end),edges)./data_split;
line(edges,bins_b,'color','k','LineWidth',3)
line(edges,bins_i,'color','r','LineWidth',3)
ylabel('p(FF)','FontSize',16)
kstest2(Result_pre.FanoFactor(1:data_split),Result_pre.FanoFactor(data_split + 1 :end))

subplot(4,1,4);hold on
FF = [(Result_pre.FanoFactor + Result_pre.Fano_95CIs(:,1)); ...
    (Result_pre.FanoFactor + Result_pre.Fano_95CIs(:,2))];
range = max(FF)-min(FF);
bins = (range / 10);
edges = [min(FF):bins:max(FF)];
FF_b = [(Result_pre.FanoFactor(1:data_split,1) + Result_pre.Fano_95CIs(1:data_split,1)); ...
    (Result_pre.FanoFactor(1:data_split,1) + Result_pre.Fano_95CIs(1:data_split,2))];
bins_b = histc(FF_b,edges)./size(FF_b,1);
FF_i = [(Result_pre.FanoFactor(data_split + 1 :end,1) + Result_pre.Fano_95CIs(data_split + 1 :end,1)); ...
    (Result_pre.FanoFactor(data_split + 1 :end,1) + Result_pre.Fano_95CIs(data_split + 1 :end,2))];
bins_i = histc(FF_i,edges)./size(FF_i,1);
line(edges,bins_b,'color','k','LineWidth',3)
line(edges,bins_i,'color','r','LineWidth',3)
ylabel('p(FF) : 95%CI added','FontSize',16)

FF_b = Result_pre.FanoFactor(1:data_split);
FF_i = Result_pre.FanoFactor(data_split + 1 :end);
[h,p] = kstest2(FF_b,FF_i)
intact_prestim_FF = median(FF_b)
CI = getCDFconf(FF_b,95)
gaba_prestim_FF = median(FF_i)
CI = getCDFconf(FF_b,95)

set(hfig,'Position',[440    92   487   706]);

%%%%%%%%%%%%%%%
%for intact condition...

% times  = 25:50:size(tmp_post,2)-24; % from 1950 ms before target onset until 2000 ms after.
% fanoParams.alignTime = (size(tmp_post,2)/2)-25; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
% fanoParams.matchReps = 10;
% fanoParams.binSpacing = 0.25;
% fanoParams.boxWidth = 50;
Result_post = VarVsMean(post_data, times, fanoParams); 
plotFano(Result_post);

hfig = figure;hold on
subplot(4,1,1)
hold on
data_split = size(Result_post.meanRateAll,1)/2;
% plot(Result_post.meanRateAll(1:data_split),'--.k','LineWidth',2)
% plot(Result_post.meanRateAll(data_split + 1 :end),'--.r','LineWidth',2)
plot(Result_post.meanRateSelect(1:data_split),'color','k','LineWidth',3)
plot(Result_post.meanRateSelect(data_split + 1 :end),'color','r','LineWidth',3)
axis tight
ylabel('mean count Raw','FontSize',16)

subplot(4,1,2);hold on
% plot(Result_post.FanoFactorAll(1:data_split),'--.k','LineWidth',2)
% plot(Result_post.FanoFactorAll(data_split + 1 :end),'--.r','LineWidth',2)
plot(Result_post.FanoFactor(1:data_split),'color','k','LineWidth',3)
plot(Result_post.FanoFactor(data_split + 1 :end),'color','r','LineWidth',3)
axis tight
ylabel('FanoFactor matched','FontSize',16)

subplot(4,1,3);hold on
range = max(Result_post.FanoFactor)-min(Result_post.FanoFactor);
bins = (range / 10);
edges = [min(Result_post.FanoFactor):bins:max(Result_post.FanoFactor)];
bins_b = histc(Result_post.FanoFactor(1:data_split),edges)./data_split;
bins_i = histc(Result_post.FanoFactor(data_split + 1 :end),edges)./data_split;
line(edges,bins_b,'color','k','LineWidth',3)
line(edges,bins_i,'color','r','LineWidth',3)
ylabel('p(FF)','FontSize',16)
kstest2(Result_post.FanoFactor(1:data_split),Result_post.FanoFactor(data_split + 1 :end))

subplot(4,1,4);hold on
FF = [(Result_post.FanoFactor + Result_post.Fano_95CIs(:,1)); ...
    (Result_post.FanoFactor + Result_post.Fano_95CIs(:,2))];
range = max(FF)-min(FF);
bins = (range / 10);
edges = [min(FF):bins:max(FF)];
FF_b = [(Result_post.FanoFactor(1:data_split,1) + Result_post.Fano_95CIs(1:data_split,1)); ...
    (Result_post.FanoFactor(1:data_split,1) + Result_post.Fano_95CIs(1:data_split,2))];
bins_b = histc(FF_b,edges)./size(FF_b,1);
FF_i = [(Result_post.FanoFactor(data_split + 1 :end,1) + Result_post.Fano_95CIs(data_split + 1 :end,1)); ...
    (Result_post.FanoFactor(data_split + 1 :end,1) + Result_post.Fano_95CIs(data_split + 1 :end,2))];
bins_i = histc(FF_i,edges)./size(FF_i,1);
line(edges,bins_b,'color','k','LineWidth',3)
line(edges,bins_i,'color','r','LineWidth',3)
ylabel('p(FF) : 95%CI added','FontSize',16)
set(hfig,'Position',[440    92   487   706]);


FF_b = Result_post.FanoFactor(1:data_split);
FF_i = Result_post.FanoFactor(data_split + 1 :end);
[h,p] = kstest2(FF_b,FF_i)
intact_stimulus_FF = median(FF_b)
CI = getCDFconf(FF_b,95)
gaba_stimulus_FF = median(FF_i)
CI = getCDFconf(FF_i,95)

%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start here to plot normalized onset and offset and fitlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if highres ==1;
    boxW = 20;
    %%%%%%%%%%%%%%
    %get offset data
    zerotime = 7000;
    %used aligned spikes trains so that end is at beginning
    for icond = 1:size(base_data_All,2)
        condlen = size(base_data_All(icond).spikes,2);
        condend = condlen - zerotime + 2000 - 1;
        lastbit(icond).spikes = base_data_All(icond).spikes(:,(condend):end);
        lastbit_gz(icond).spikes = gz_data_All(icond).spikes(:,(condend):end);
    end
    times = (zerotime-4000)+boxW/2:boxW/2:(8999-boxW); % from 4000 ms before target offset until 2000 ms after.
    fanoParams.alignTime = zerotime-(boxW/2); % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
    fanoParams.matchReps = 10;
    fanoParams.binSpacing = 0.25;
    fanoParams.boxWidth = boxW;
    Result_base_end = VarVsMean(lastbit, times, fanoParams);
%     plotFano(Result_base_end);
    
    Result_gz_end = VarVsMean(lastbit_gz, times, fanoParams);
    % plotFano(Result_gz_end);
    
    figure;hold on
    line(Result_base_end.times,Result_base_end.meanRateAll,'color','k')
    line(Result_gz_end.times,Result_gz_end.meanRateAll,'color','r')
    set(gca,'XLim',[-1000,1500])
    figure;hold on
    line(Result_base_end.times,Result_base_end.FanoFactorAll,'color','k')
    line(Result_gz_end.times,Result_gz_end.FanoFactorAll,'color','r')
    
    %%%%%%%%%%%%%%%%
    %get onset data
    zerotime = 2000;
    times = boxW/2:boxW/2:zerotime + 4000 - boxW/2; % from 1950 ms before target onset until 2000 ms after.
    fanoParams.alignTime = zerotime - (boxW/2); % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
    fanoParams.matchReps = 10;
    fanoParams.binSpacing = 0.25;
    fanoParams.boxWidth = boxW;
    Result_base = VarVsMean(base_data_All, times, fanoParams);
    % plotFano(Result_base);
    
    Result_gz = VarVsMean(gz_data_All, times, fanoParams);
    % plotFano(Result_gz);
    figure;hold on
    line(Result_base.times,Result_base.meanRateAll,'color','k')
    line(Result_gz.times,Result_gz.meanRateAll,'color','r')
    axis tight
    figure;hold on
    plot(Result_base.times,Result_base.FanoFactor,'color','k')
    plot(Result_base.times,Result_gz.FanoFactor,'color','r')
    axis tight
end

if highres ==0
    %%%%%%%%%%%%%%
    %get offset data
    
    %used aligned spikes trains so that end is at beginning
    for icond = 1:size(base_data_All,2)
        condlen = size(base_data_All(icond).spikes,2);
        condend = condlen - 8999;
        lastbit(icond).spikes = base_data_All(icond).spikes(:,(condend):end);
        lastbit_gz(icond).spikes = gz_data_All(icond).spikes(:,(condend):end);
    end
    times = 3025:50:8950; % from 5000 ms before target onset until 2000 ms after.
    fanoParams.alignTime = 6975; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
    fanoParams.matchReps = 10;
    fanoParams.binSpacing = 0.25;
    fanoParams.boxWidth = 50;
    Result_base_end = VarVsMean(lastbit, times, fanoParams);
    % plotFano(Result_base_end);
    
    Result_gz_end = VarVsMean(lastbit_gz, times, fanoParams);
    % plotFano(Result_gz_end);
    
    figure;hold on
    line(Result_base_end.times,Result_base_end.meanRateAll,'color','k')
    line(Result_gz_end.times,Result_gz_end.meanRateAll,'color','r')
    set(gca,'XLim',[-1000,1500])
    figure;hold on
    line(Result_base_end.times,Result_base_end.FanoFactorAll,'color','k')
    line(Result_gz_end.times,Result_gz_end.FanoFactorAll,'color','r')
    
    %%%%%%%%%%%%%%%%
    %get onset data
    
    times = 25:50:6000; % from 1950 ms before target onset until 2000 ms after.
    fanoParams.alignTime = 1975; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
    fanoParams.matchReps = 10;
    fanoParams.binSpacing = 0.25;
    fanoParams.boxWidth = 50;
    Result_base = VarVsMean(base_data_All, times, fanoParams);
    % plotFano(Result_base);
    
    Result_gz = VarVsMean(gz_data_All, times, fanoParams);
    % plotFano(Result_gz);
    figure;hold on
    line(Result_base.times,Result_base.meanRateAll,'color','k')
    line(Result_gz.times,Result_gz.meanRateAll,'color','r')
    axis tight
    figure;hold on
    plot(Result_base.times,Result_base.FanoFactor,'color','k')
    plot(Result_base.times,Result_gz.FanoFactor,'color','r')
    axis tight
    
end


%%%%%%%%%%%%%%%
% plot onset and offset as function of pre-stim FF

% onset as a % of mean baseline FF and quantify
onsetind = max(find(Result_base.times <= 0));
prestimFF_base = mean(Result_base.FanoFactor(1:onsetind));
prestimFF_gz = mean(Result_gz.FanoFactor(1:onsetind));
norm_base = Result_base.FanoFactor / prestimFF_base;
norm_base_CI = Result_base.Fano_95CIs / prestimFF_base;
meanDropFF_base = mean(norm_base(onsetind+5:end));
norm_gaba = Result_gz.FanoFactor / prestimFF_gz;
norm_gaba_CI = Result_gz.Fano_95CIs/ prestimFF_gz;
meanDropFF_gaba = mean(norm_gaba(onsetind+5:end));

%  offset  as a % pre-stim mean baseline FF
offsetind = max(find(Result_base_end.times < 0));
norm_base_end = Result_base_end.FanoFactor / prestimFF_base;
norm_base_CI_end = Result_base_end.Fano_95CIs / prestimFF_base;
meanDropFFend_base = mean(norm_base_end(5:offsetind-5));
norm_gaba_end = Result_gz_end.FanoFactor / prestimFF_gz;
norm_gaba_CI_end = Result_gz_end.Fano_95CIs / prestimFF_gz;
meanDropFFend_gaba = mean(norm_gaba_end(5:offsetind-5));

avgdrop_base = mean([meanDropFF_base,meanDropFFend_base])
avgdrop_gz = mean([meanDropFF_gaba,meanDropFFend_gaba])

% Fano in intact condition aligned at on and off
%also get mean Fano relative to baseline and CI
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000);
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000);
line(onsetT(1:clipend),norm_base(1:clipend),'color','k','LineWidth',2)
line(offsetT(clipon:end),norm_base_end(clipon:end),'color','k','LineWidth',2)
% line(onsetT(1:clipend),norm_gaba(1:clipend),'color','r')
% line(offsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
line([0,0],[0.6,1],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0.6,1],'color','r')
line(onsetT(1:clipend), norm_base_CI(1:clipend,1),'color',[0.5 0.5 0.5])
line(onsetT(1:clipend), norm_base_CI(1:clipend,2),'color',[0.5 0.5 0.5])
line(offsetT(clipon:end),norm_base_CI_end(clipon:end,1),'color',[0.5 0.5 0.5])
line(offsetT(clipon:end),norm_base_CI_end(clipon:end,2),'color',[0.5 0.5 0.5])

medianFano = median(norm_base(onsetind:clipend))
CI = getCDFconf(norm_base(onsetind:clipend),95)
medianCIacrossTime = median(norm_base_CI(onsetind:clipend,:))
numconditions = size(Result_pre.meanRateAll_Trials,1)
[h,p]=kstest([norm_base(onsetind:clipend)',norm_base(1:onsetind)'])
[h,p]=kstest([Result_base.FanoFactor(onsetind:clipend)',Result_base.FanoFactor(1:onsetind)'])


% make figure for onset align and offset align of mean and FF using % baseline not absolute
%plotted on same axes for making figure
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000);
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000);
line(onsetT(1:clipend),norm_base(1:clipend),'color','k')
line(offsetT(clipon:end),norm_base_end(clipon:end),'color','k')
line(onsetT(1:clipend),norm_gaba(1:clipend),'color','r')
line(offsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
line([0,0],[0.6,1],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0.6,1],'color','r')
% saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormFF.fig')

%mean rates across both conditions plotted on same axes for making figure
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000);
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000);
line(onsetT(1:clipend),Result_base.meanRateAll(1:clipend))
line(offsetT(clipon:end),Result_base_end.meanRateAll(clipon:end))
line(onsetT(1:clipend),Result_gz.meanRateAll(1:clipend),'color','r')
line(offsetT(clipon:end),Result_gz_end.meanRateAll(clipon:end),'color','r')
line([0,0],[0,12],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0,12],'color','r')
prestimRate = mean(Result_base.meanRateAll(1:onsetind))
% saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormRates.fig')

%normalized to max for both conditions at stim onset and offset 
%for figure on same time course for intact and gabazine
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000);
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000);
line(onsetT(1:clipend),Result_base.FanoFactor(1:clipend)/max(Result_base.FanoFactor(1:clipend)),'color','k')
line(offsetT(clipon:end),Result_base_end.FanoFactor(clipon:end)/max(Result_base_end.FanoFactor(clipon:end)),'color','k')
line(onsetT(1:clipend),Result_gz.FanoFactor(1:clipend)/max(Result_gz.FanoFactor(1:clipend)),'color','r')
line(offsetT(clipon:end),Result_gz_end.FanoFactor(clipon:end)/max(Result_gz_end.FanoFactor(clipon:end)),'color','r')
line([0,0],[0.6,1],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0.6,1],'color','r')

%scatter of mean rates to get a a sense of scale and offset
figure;hold on
scatter(Result_base.meanRateAll,Result_gz.meanRateAll)
scatter(Result_base_end.meanRateAll,Result_gz_end.meanRateAll)
line([0,25],[0,25],'color','r')
axis tight
axis image

%%%%%%%%%%%%
% get fits for onset and offset under both conditions

%onset
options = fitoptions('exp1');
lowerTau = -1/20;
upperTau = -1/700;
options.Lower = [0,lowerTau];
options.Upper = [Inf,upperTau];
fitT = Result_base.times(onsetind:end);
fitY = norm_base(onsetind:end) - meanDropFF_base;
figure;hold on
line(fitT,fitY)
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
line(fitT,fitline,'color','r')
onsetTau = -1/cf.b
figure;hold on
line(Result_base.times,norm_base,'color','k')
line(fitT,fitline + meanDropFF_base,'color','r')
% set(gca,'XLim',[-500,500],'YLim',[0.6,1.2])

%offset
options = fitoptions('exp1');
lowerTau = 1/20;
upperTau = 1/700;
options.Lower = [0,-Inf];
options.Upper = [Inf,Inf];
stopfit = find(Result_base_end.times == 500)
fitT = Result_base_end.times(1:stopfit);
fitY = norm_base_end(1:stopfit) - meanDropFFend_base;
figure;hold on
line(fitT,fitY)
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
line(fitT,fitline,'color','r')
offsetTau = -1/cf.b
figure;hold on
line(Result_base_end.times,norm_base_end,'color','k')
line(fitT,fitline + meanDropFFend_base,'color','r')
% set(gca,'XLim',[-500,1000],'YLim',[0.5,1.5])

% make figure for onset align and offset align of mean and FF using % baseline not absolute
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000);
line(onsetT(1:clipend),norm_base(1:clipend),'color','k')
line(onsetT(1:clipend),norm_gaba(1:clipend),'color','r')
%onset
options = fitoptions('exp1');
lowerTau = -1/20;
upperTau = -1/700;
options.Lower = [0,lowerTau];
options.Upper = [Inf,upperTau];
fitT = Result_base.times(onsetind:end);
fitY = norm_base(onsetind:end) - meanDropFF_base;
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
onsetTau = -1/cf.b
line(fitT,fitline + meanDropFF_base,'color','g')

offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000);
line(offsetT(clipon:end),norm_base_end(clipon:end),'color','k')
line(offsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
%offset
options = fitoptions('exp1');
options.Lower = [0,-Inf];
options.Upper = [Inf,Inf];
stopfit = find(Result_base_end.times == 500)
fitT = Result_base_end.times(1:stopfit);
fitY = norm_base_end(1:stopfit) - meanDropFFend_base;
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
offsetTau = -1/cf.b
line(fitT+size(Result_base.times,1)*fanoParams.boxWidth,fitline + meanDropFFend_base,'color','g')

line([0,0],[0.6,1],'color','b')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0.6,1],'color','b')
axis tight
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormFF_on_off_withFIT.fig')
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormFF_on_off_withFIT.jpg')


%onset gaba
options = fitoptions('exp1');
lowerTau = -1/20;
upperTau = -1/700;
options.Lower = [0,lowerTau];
options.Upper = [Inf,upperTau];
fitT = Result_gz.times(onsetind:end);
fitY = norm_gaba(onsetind:end) - meanDropFF_gaba;
figure;hold on
line(fitT,fitY,'color','r')
cf = fit(fitT,fitY,'exp1',options);
fitline_on_gaba = cf.a*exp(cf.b*fitT);
line(fitT,fitline_on_gaba,'color','g')
onsetTau_gaba = -1/cf.b

%offset gaba
options = fitoptions('exp1');
options.Lower = [0,-Inf];
options.Upper = [Inf,Inf];
stopfit = find(Result_gz_end.times == 500)
fitT_end_gaba = Result_gz_end.times(1:stopfit);
fitY = norm_gaba_end(1:stopfit) - meanDropFFend_gaba;
figure;hold on
line(fitT_end_gaba,fitY,'color','r')
cf = fit(fitT_end_gaba,fitY,'exp1',options);
fitline_off_gaba = cf.a*exp(cf.b*fitT_end_gaba);
line(fitT_end_gaba,fitline_off_gaba,'color','g')
offsetTau_gaba = -1/cf.b
% 
% 
% hfig = figure;
% hold on
% onsetT = Result_base_end.times;
% clipon = find(Result_base_end.times == -2000);
% line(onsetT(clipon:end),norm_base_end(clipon:end),'color','k')
% line(onsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
% %onset
% options = fitoptions('exp1');
% lowerTau = -1/20;
% upperTau = -1/700;
% options.Lower = [0,lowerTau];
% options.Upper = [Inf,upperTau];
% fitT = Result_gz.times(onsetind:end);
% fitY = norm_gaba(onsetind:end) - meanDropFF_gaba;
% cf = fit(fitT,fitY,'exp1',options);
% fitline_on_gaba = cf.a*exp(cf.b*fitT);
% onsetTau_gaba = -1/cf.b
% line(fitT,fitline_on_gaba + meanDropFF_gaba,'color','g')
% 
% offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
% clipon = find(Result_base_end.times == -2000);
% line(offsetT(clipon:end),norm_base_end(clipon:end),'color','k')
% line(offsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
% %offset
% options = fitoptions('exp1');
% options.Lower = [0,-Inf];
% options.Upper = [Inf,Inf];
% stopfit = find(Result_base_end.times == 700)
% fitT = Result_base_end.times(1:stopfit);
% fitY = norm_base_end(1:stopfit) - meanDropFFend_base;
% cf = fit(fitT,fitY,'exp1',options);
% fitline = cf.a*exp(cf.b*fitT);
% offsetTau = -1/cf.b
% line(fitT+size(Result_base.times,1)*fanoParams.boxWidth,fitline + meanDropFFend_base,'color','g')
% 
% line([0,0],[0.6,1],'color','b')
% line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
%     [0.6,1],'color','b')
% axis tight

%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%% Do the FanoFactor analysis across all cells...

times = 25:50:6000; % from 1950 ms before target onset until 2000 ms after.
fanoParams.alignTime = 1975; % this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
fanoParams.matchReps = 10;
fanoParams.binSpacing = 0.25;
fanoParams.boxWidth = 50;
Result_base = VarVsMean(base_data_All, times, fanoParams); 
plotFano(Result_base);

figure;hold on
line(Result_base.times,Result_base.FanoFactor','color','r')
line(Result_base.times,Result_base.FanoFactorAll','color','k')
legend('meanMatched','unmatched')
ylabel('FanoFactor', 'FontSize',14);
xlabel('milliseconds','FontSize',14);
% mv = ScatterMovie(Result_base);
% yourName = '/Users/kperks/GitHub/iontoncm/FANOanalysis/baseline_condition_matched.avi';
% movie2avi(mv, yourName, 'FPS', 15, 'compression', 'none');

Result_gz = VarVsMean(gz_data_All, times, fanoParams); 
plotFano(Result_gz);

figure;hold on
line(Result_gz.times,Result_gz.FanoFactor','color','r')
line(Result_gz.times,Result_gz.FanoFactorAll','color','k')
legend('meanMatched','unmatched')
ylabel('FanoFactor', 'FontSize',14);
xlabel('milliseconds','FontSize',14);
% mv = ScatterMovie(Result_base);
% yourName = '/Users/kperks/GitHub/iontoncm/FANOanalysis/InhibitionBlocked_condition_matched.avi';
% movie2avi(mv, yourName, 'FPS', 15, 'compression', 'none');


figure;hold on
line(Result_base.times,Result_base.meanRateAll,'color','k')
line(Result_gz.times,Result_gz.meanRateAll,'color','r')
axis tight

% 
% figure;hold on
% plot(Result_base.FanoFactorAll,'color','k')
% plot(Result_gz.FanoFactorAll,'color','r')
% axis tight
% 
% %when did mean matching...
% figure;hold on
% plot(Result_base.meanRateSelect,'color','k')
% plot(Result_gz.meanRateSelect,'color','r')
% axis tight


figure;hold on
plot(Result_base.times,Result_base.FanoFactor,'color','k')
plot(Result_base.times,Result_gz.FanoFactor,'color','r')
axis tight

% plot as a % of mean baseline FF and quantify
onsetind = max(find(Result_base.times <= 0))
prestimFF_base = mean(Result_base.FanoFactor(1:onsetind));
prestimFF_gz = mean(Result_gz.FanoFactor(1:onsetind));
norm_base = Result_base.FanoFactor / mean(Result_base.FanoFactor(1:onsetind));
meanDropFF_base = mean(norm_base(onsetind+5:end));
norm_gaba = Result_gz.FanoFactor / mean(Result_gz.FanoFactor(1:onsetind));
meanDropFF_gaba = mean(norm_gaba(onsetind+5:end));
figure;hold on
plot(Result_base.times , norm_base,'color','k')
plot(Result_base.times , norm_gaba,'color','r')

norm_base_onsetAlign = norm_base;


% plot aligned to offset results as % pre-stim mean
offsetind = max(find(Result_base_end.times < 0))
norm_base_end = Result_base_end.FanoFactor / prestimFF_base;
meanDropFFend_base = mean(norm_base_end(5:offsetind-5));
norm_gaba_end = Result_gz_end.FanoFactor / prestimFF_gz;
meanDropFFend_gaba = mean(norm_gaba_end(5:offsetind-5));
figure;hold on
plot(Result_base_end.times , norm_base_end,'color','k')
plot(Result_base_end.times , norm_gaba_end,'color','r')

avgdrop_base = mean([meanDropFF_base,meanDropFFend_base])
avgdrop_gz = mean([meanDropFF_gaba,meanDropFFend_gaba])
%% re-make figure for onset align and offset align of mean and FF using % baseline not absolute
hfig = figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000)
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000)
line(onsetT(1:clipend),norm_base(1:clipend),'color','k')
line(offsetT(clipon:end),norm_base_end(clipon:end),'color','k')
line(onsetT(1:clipend),norm_gaba(1:clipend),'color','r')
line(offsetT(clipon:end),norm_gaba_end(clipon:end),'color','r')
line([0,0],[0.6,1],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0.6,1],'color','r')
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormFF.fig')

figure;
hold on
onsetT = Result_base.times;
clipend = find(Result_base.times == 2000)
offsetT = Result_base_end.times + size(Result_base.times,1)*fanoParams.boxWidth;
clipon = find(Result_base_end.times == -2000)
line(onsetT(1:clipend),Result_base.meanRateAll(1:clipend))
line(offsetT(clipon:end),Result_base_end.meanRateAll(clipon:end))
line(onsetT(1:clipend),Result_gz.meanRateAll(1:clipend),'color','r')
line(offsetT(clipon:end),Result_gz_end.meanRateAll(clipon:end),'color','r')

line([0,0],[0,12],'color','r')
line([size(Result_base.times,1)*fanoParams.boxWidth,size(Result_base.times,1)*fanoParams.boxWidth],...
    [0,12],'color','r')
saveas(hfig,'/Volumes/kperks/Documents/EJN/Figures/EPSforFIG4/NormRates.fig')


figure;hold on
scatter(Result_base.meanRateAll,Result_gz.meanRateAll)
scatter(Result_base_end.meanRateAll,Result_gz_end.meanRateAll)
line([0,25],[0,25],'color','r')
axis tight
axis image

prestimRate = mean(Result_base.meanRateAll(1:onsetind))

%% fit normalized (by baseline mean) FF curves to get tau

%onset
options = fitoptions('exp1');
lowerTau = -1/20;
upperTau = -1/700;
options.Lower = [0,lowerTau];
options.Upper = [Inf,upperTau];
fitT = Result_base.times(onsetind:end);
fitY = norm_base(onsetind:end) - meanDropFF_base;
figure;hold on
line(fitT,fitY)
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
line(fitT,fitline,'color','r')
onsetTau = -1/cf.b
figure;hold on
line(Result_base.times,norm_base,'color','k')
line(fitT,fitline + meanDropFF_base,'color','r')
set(gca,'XLim',[-500,500],'YLim',[0.6,1.2])

%offset
options = fitoptions('exp1');
lowerTau = 1/20;
upperTau = 1/700;
options.Lower = [0,-Inf];
options.Upper = [Inf,Inf];
stopfit = find(Result_base_end.times == 700)
fitT = Result_base_end.times(1:stopfit);
fitY = norm_base_end(1:stopfit) - meanDropFFend_base;
figure;hold on
line(fitT,fitY)
cf = fit(fitT,fitY,'exp1',options);
fitline = cf.a*exp(cf.b*fitT);
line(fitT,fitline,'color','r')
offsetTau = -1/cf.b
set(gca,'XLim',[-500,1000],'YLim',[0.5,1.5])

%onset gaba
options = fitoptions('exp1');
lowerTau = -1/20;
upperTau = -1/700;
options.Lower = [0,lowerTau];
options.Upper = [Inf,upperTau];
fitT = Result_gz.times(onsetind:end);
fitY = norm_gaba(onsetind:end) - meanDropFF_gaba;
figure;hold on
line(fitT,fitY,'color','r')
cf = fit(fitT,fitY,'exp1',options);
fitline_on_gaba = cf.a*exp(cf.b*fitT);
line(fitT,fitline,'color','g')
onsetTau_gaba = -1/cf.b

%offset gaba
options = fitoptions('exp1');
options.Lower = [0,-Inf];
options.Upper = [Inf,Inf];
stopfit = find(Result_gz_end.times == 450)
fitT_end_gaba = Result_gz_end.times(1:stopfit);
fitY = norm_gaba_end(1:stopfit) - meanDropFFend_gaba;
figure;hold on
line(fitT_end_gaba,fitY,'color','r')
cf = fit(fitT_end_gaba,fitY,'exp1',options);
fitline_off_gaba = cf.a*exp(cf.b*fitT_end_gaba);
line(fitT_end_gaba,fitline_off_gaba,'color','g')
offsetTau_gaba = -1/cf.b
%%
zeroind = 0 %(because bin inds are the center of each bin, 
...the last bin with only pre-stim data is the one half a bin size below zero
%
% find time constant for both onset and offset 
% find best fit by increasing time included and finding max r-sq

% stimulus offset... work forwards in samps


%****Result_base_end
xtime = Result_base_end.times;
target = find(xtime == zeroind);
fano_offset =  median(Result_base_end.FanoFactorAll(1:target,:));
ffend = Result_base_end.FanoFactorAll - fano_offset ;
sampind = 1;
rsq = [];
for isamp = target:size(ffend,1)
ffend_clip = ffend(1:isamp,1);
xclip = xtime(1:isamp,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
[rsq(isamp) rmse] = rsquare(ffend_clip,fitline);
sampind = sampind + 1;
end
maxfitind = find(rsq == max(rsq));
ffend_clip = ffend(1:maxfitind,1);
xclip = xtime(1:maxfitind,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
tao_end = -1/cf.b;

hfig = figure;
hold on
line(xtime,ffend + fano_offset,'color','k','LineWidth',3)
line(xclip,fitline + fano_offset,'color','r','LineWidth',2)
axis tight
ylims = get(gca,'YLim');
text(-4000,mean(ylims),['Time constant (msec) = ' num2str(round(tao_end))]...
    ,'Interpreter','none','FontSize',14)
xlabel('milliseconds','FontSize',14)
ylabel('Fano Factor','FontSize',14)
line([0,0],ylims,'color','k')
text(0,ylims(2),'stim offset','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%****Result_gz_end
xtime = Result_gz_end.times;
target = find(xtime == zeroind);
fano_offset =  median(Result_gz_end.FanoFactorAll(1:target,:));
ffend = Result_gz_end.FanoFactorAll - fano_offset ;
sampind = 1;
rsq = [];
for isamp = target:size(ffend,1)
ffend_clip = ffend(1:isamp,1);
xclip = xtime(1:isamp,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
[rsq(isamp) rmse] = rsquare(ffend_clip,fitline);
sampind = sampind + 1;
end
maxfitind = find(rsq == max(rsq));
ffend_clip = ffend(1:maxfitind,1);
xclip = xtime(1:maxfitind,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
tao_end = -1/cf.b;

hfig = figure;
hold on
line(xtime,ffend + fano_offset,'color','k','LineWidth',3)
line(xclip,fitline + fano_offset,'color','r','LineWidth',2)
axis tight
ylims = get(gca,'YLim');
text(-4000,mean(ylims),['Gabazine Time constant (msec) = ' num2str(round(tao_end))]...
    ,'Interpreter','none','FontSize',14)
xlabel('milliseconds','FontSize',14)
ylabel('Fano Factor','FontSize',14)
line([0,0],ylims,'color','k')
text(0,ylims(2),'stim offset','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%plot both FANO traces normalized together
hfig = figure;
hold on
xtime = Result_gz_end.times;
target = find(xtime == zeroind);
fano_offset =  median(Result_gz_end.FanoFactorAll(1:target,:));
norm_gz_end = Result_gz_end.FanoFactorAll - fano_offset;
norm_gz_end = norm_gz_end / max(norm_gz_end);
line(xtime,norm_gz_end,'color','r','LineWidth',3)
xtime = Result_base_end.times;
target = find(xtime == zeroind);
fano_offset =  median(Result_base_end.FanoFactorAll(1:target,:));
norm_base_end = Result_base_end.FanoFactorAll - fano_offset;
norm_base_end = norm_base_end / max(norm_base_end);
line(xtime,norm_base_end,'color','b','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
legend('gabazine','baseline Pharm')
line([0,0],ylims,'color','k')
xlabel('milliseconds','FontSize',14)
ylabel('normalized FanoFactor','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%****fit normalized gz data and base data
%************did not change Tao estimate


options = fitoptions('exp1');
% stim onset... flip and then work forward in samps.. then flip back to
% baseline Pharm condition
xtime = Result_base.times;
ffend = Result_base.FanoFactorAll;
ffend_flip = fliplr(ffend')';
xtime_flip = -fliplr(xtime')';
target = find(xtime_flip == zeroind);
fano_offset =  median(ffend_flip(1:target,:));
ffend_flip = ffend_flip - fano_offset;

rsq = [];
for isamp = target-10:target+5 ; %size(ffend_flip,1)
ffend_clip = ffend_flip(1:isamp,1);
xclip = xtime_flip(1:isamp,1);
options.Lower = [0,0];
options.Upper = [Inf,Inf];
cf = fit(xclip,ffend_clip,'exp1',options);
fitline = cf.a*exp(cf.b*xclip);
[rsq(isamp) rmse] = rsquare(ffend_clip,fitline);
end
maxfitind = find(rsq == max(rsq));
ffend_clip = ffend_flip(1:maxfitind,1);
xclip = xtime_flip(1:maxfitind,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
tao_end = -1/cf.b;

hfig = figure;
hold on
line(-fliplr(xtime_flip')',fliplr(ffend_flip')' + fano_offset,'color','k','LineWidth',3)
line(-fliplr(xclip')',fliplr((fitline + fano_offset)')','color','r','LineWidth',2)
axis tight
ylims = get(gca,'YLim');
text(2000,mean(ylims),['Time constant (msec) = ' num2str(round(tao_end))]...
    ,'Interpreter','none','FontSize',14)
xlabel('milliseconds','FontSize',14)
ylabel('Fano Factor','FontSize',14)
line([0,0],ylims,'color','k')
text(0,ylims(2),'stim onset','FontSize',14)
set(hfig,'Position',[388         517    1038         281])


% gabazine Pharm condition
xtime = Result_gz.times;
ffend = Result_gz.FanoFactorAll;
ffend_flip = fliplr(ffend')';
xtime_flip = -fliplr(xtime')';
target = find(xtime_flip == 0);
fano_offset =  median(ffend_flip(1:target,:));
ffend_flip = ffend_flip - fano_offset;

rsq = [];
for isamp = target-10:target+5 ; %size(ffend_flip,1)
ffend_clip = ffend_flip(1:isamp,1);
xclip = xtime_flip(1:isamp,1);
options.Lower = [0,0];
options.Upper = [Inf,Inf];
cf = fit(xclip,ffend_clip,'exp1',options);
fitline = cf.a*exp(cf.b*xclip);
[rsq(isamp) rmse] = rsquare(ffend_clip,fitline);
end
maxfitind = find(rsq == max(rsq));
ffend_clip = ffend_flip(1:maxfitind,1);
xclip = xtime_flip(1:maxfitind,1);
cf = fit(xclip,ffend_clip,'exp1');
fitline = cf.a*exp(cf.b*xclip);
tao_end = -1/cf.b;

hfig = figure;
hold on
line(-fliplr(xtime_flip')',fliplr(ffend_flip')' + fano_offset,'color','k','LineWidth',3)
line(-fliplr(xclip')',fliplr((fitline + fano_offset)')','color','r','LineWidth',2)
axis tight
ylims = get(gca,'YLim');
text(2000,mean(ylims),['Gabazine Time constant (msec) = ' num2str(round(tao_end))]...
    ,'Interpreter','none','FontSize',14)
xlabel('milliseconds','FontSize',14)
ylabel('Fano Factor','FontSize',14)
line([0,0],ylims,'color','k')
text(0,ylims(2),'stim onset','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%plot both FANO traces normalized together
hfig = figure;
hold on
xtime = Result_gz.times;
target = find(xtime == 0);
fano_offset =  median(Result_gz.FanoFactorAll(target:end,:));
norm_gz = Result_gz.FanoFactorAll - fano_offset;
norm_gz = norm_gz / max(norm_gz);
line(xtime,norm_gz,'color','r','LineWidth',3)
xtime = Result_base.times;
target = find(xtime == 0);
fano_offset =  median(Result_base.FanoFactorAll(target:end,:));
norm_base = Result_base.FanoFactorAll - fano_offset;
norm_base = norm_base / max(norm_base);
line(xtime,norm_base,'color','b','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
legend('gabazine','baseline Pharm')
line([0,0],ylims,'color','k')
xlabel('milliseconds','FontSize',14)
ylabel('normalized FanoFactor','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%% Spike Rates "GAIN" intact and gabazine...
%%%%%%%%%% test if multiplicative gain...
%responses aligned to stim onset
figure;
scatter(Result_base.meanRateAll,Result_gz.meanRateAll,30,'k')
options = fitoptions('poly1');
options.Upper = [Inf,0];
options.Lower = [-Inf,0];
cf = fit(Result_base.meanRateAll,Result_gz.meanRateAll,'poly1',options)
fitline = cf.p1*Result_base.meanRateAll;
Rsq = rsquare(Result_gz.meanRateAll,fitline);
line(Result_base.meanRateAll,fitline,'color','r')
set(gca,'YLim',[0,25],'XLim',[0,25])
text(5,1,['gain = ' num2str(cf.p1) '  R^2 = ' num2str(Rsq)],'FontSize',14)
xlabel('spike rate baseline','FontSize',14)
ylabel('spike rate gabazine','FontSize',14)
title('trials aligned to stimulus onset','FontSize',14)
minPreStim = min(Result_base.meanRateAll(1:find(Result_base.times == 0)));
maxPreStim = max(Result_base.meanRateAll(1:find(Result_base.times == 0)));
line([minPreStim,minPreStim],[0,25])
line([maxPreStim,maxPreStim],[0,25])


%responses aligned to stim offset
figure;
scatter(Result_base_end.meanRateAll,Result_gz_end.meanRateAll,30,'k')
options = fitoptions('poly1');
options.Upper = [Inf,0];
options.Lower = [-Inf,0];
cf = fit(Result_base_end.meanRateAll,Result_gz_end.meanRateAll,'poly1',options)
fitline = cf.p1*Result_base_end.meanRateAll;
Rsq = rsquare(Result_gz_end.meanRateAll,fitline);
line(Result_base_end.meanRateAll,fitline,'color','r')
set(gca,'YLim',[0,25],'XLim',[0,25])
text(5,1,['gain = ' num2str(cf.p1) '  R^2 = ' num2str(Rsq)],'FontSize',14)
xlabel('spike rate baseline','FontSize',14)
ylabel('spike rate gabazine','FontSize',14)
title('trials aligned to stimulus offset','FontSize',14)

%do for each cell individually and split up pre-stim from stim
% need to do it off of the psth result, not the fano result?
% or just need to tap into the fano result
%%%%%%
%smaller boxWidth gives better linear fit for slope,
%but bigger boxWidth makes sure that changes in onset or offset of
%responses does not "falsely" impact the change in rate measure
%%%%%%%%
binwin = 100;
times = binwin/2:binwin:7000-binwin/2; % from 1950 ms before target onset until 2000 ms after.
fanoParams.alignTime = 2000-binwin/2; % (2000, which is actual start, minus half boxwidth
 %this time will become zero time fanoParams.boxWidth = 50; % 50 ms sliding window.
fanoParams.matchReps = 10;
fanoParams.binSpacing = 0.25;
fanoParams.boxWidth = 20;
Result_base = VarVsMean(base_data_All, times, fanoParams); 
Result_gz = VarVsMean(gz_data_All, times, fanoParams);
%want to add output of "mean rate all trials" to this function
plotFano(Result_base);
zeroind = crossing(Result_base.times,[],0);


figure;
hold on
line(Result_base.times,mean(Result_base.meanRateAll_Trials),'color','k')
line(Result_gz.times,mean(Result_gz.meanRateAll_Trials),'color','r')
title('mean across pop')

figure;
hold on
line(Result_base.times,median(Result_base.meanRateAll_Trials),'color','k')
line(Result_gz.times,median(Result_gz.meanRateAll_Trials),'color','r')
title('median across pop')

%try just the mean of the stim period like did for pre-stim period
%need to do for each cell because different clusters of slopes i think...
base_stim_mn = [];
gz_stim_mn = [];
p_emerge_stim = [];
p_retract_stim = [];
usebase = [];
pre_slope = [];
fitslope = [];
slopes = [];
Rsq = [];
slopeind = 1;
figure;hold on
for icond = 1:size(Result_base.meanRateAll_Trials,1)
    
    %baseline pharm
    base_trial = Result_base.meanRateAll_Trials(icond,:);
    %stim
    base_stim = double(base_trial(1,zeroind+1:end));
    base_stim_mn(icond) = mean(base_trial(1,zeroind+1:end));
    %gz pharm
    gz_trial = Result_gz.meanRateAll_Trials(icond,:);
    %stim
    gz_stim = double(gz_trial(1,zeroind+1:end));
    gz_stim_mn(icond) = mean(gz_trial(1,zeroind+1:end));
    
    
    norespinds = find(base_stim == 0);
    p_emerge_stim(icond) = size(find(gz_stim(norespinds)~=0),2) ./ size(norespinds,2);
    respinds = find(base_stim ~=0);
    p_retract_stim(icond) = size(find(gz_stim(respinds) == 0),2) / size(respinds,2);

    %pre-stim slopes
    pre_slopes(icond) =  mean(gz_trial(1,1:zeroind))/mean(base_trial(1,1:zeroind));
    baseresp = base_stim(respinds);
    gzresp = gz_stim(respinds);
    if max(size(respinds))>1
        cf = fit(baseresp',gzresp','poly1');
        fitline = cf.p1*baseresp+cf.p2;
        Rsq(icond) = rsquare(gzresp,fitline);
%         if Rsq(icond) > 0.55
            usebase(slopeind) = icond;
            slopes(slopeind) = mean(gzresp) / mean(baseresp);
            fitslope(slopeind) = cf.p1;
            slopeind = slopeind +1;
            scatter(baseresp,gzresp,30,'k')
%         end
    end
end
ylims = get(gca,'YLim');
xlims = get(gca,'XLim');
line([0,max([ylims(2),xlims(2)])],[0,max([ylims(2),xlims(2)])])
title('no response windows during baseline cond eliminated')

figure;
scatter(base_stim_mn,gz_stim_mn,30,'k','fill')
ylims = get(gca,'YLim');
xlims = get(gca,'XLim');
line([0,max([ylims(2),xlims(2)])],[0,max([ylims(2),xlims(2)])])
cf = fit(base_stim_mn',gz_stim_mn','poly1');
fitline = cf.p1*base_stim_mn+cf.p2;
rsq = rsquare(gz_stim_mn,fitline);
line(base_stim_mn,fitline,'color','r')
text(1,0.3,['gain = ' num2str(cf.p1) '  R^2 = ' num2str(rsq)],'FontSize',14)
title('all conds; mean first 5sec of stim')

figure;
hold on
edges = [0:0.1:1];
n = histc(p_retract_stim,edges);
stairs(edges,n/size(p_retract_stim,2),'color','k');
n = histc(p_emerge_stim,edges);
stairs(edges,n/size(p_emerge_stim,2),'color','r');
xlabel(['median emerge=' num2str(median(p_emerge_stim)) ...
' ; median retract =' num2str(median(p_retract_stim))])

figure;
edges = [0:0.2:10];
n = histc(slopes,edges);
bar(edges,n/size(slopes,2));
title('slopes from mean response')
xlabel(['mean=' num2str(mean(slopes)) ' ; median=' num2str(median(slopes))])

figure;
edges = [0:0.2:10];
n = histc(fitslope,edges);
bar(edges,n/size(fitslope,2));
title('slopes from fit')
xlabel(['mean=' num2str(mean(fitslope)) ' ; median=' num2str(median(fitslope))])

figure;
edges = [0:0.2:10];
n = histc(pre_slopes(usebase),edges);
bar(edges,n/size(usebase,2));
title('slopes from pre-stim with significant stim fitline')
xlabel(['mean=' num2str(mean(pre_slopes(usebase))) '; median=' num2str(median(pre_slopes(usebase)))])


base_pre = [];
base_stim = [];
gz_pre = [];
gz_stim = [];
for icond = 1:size(Result_base.meanRateAll_Trials,1)
    
    %baseline pharm
        base_trial = Result_base.meanRateAll_Trials(icond,:);
    %pre-stim
        base_pre = [base_pre, mean(base_trial(1,1:zeroind))];
    %stim
        base_stim = [base_stim, base_trial(1,zeroind+1:end)];
    %gz pharm
        gz_trial = Result_gz.meanRateAll_Trials(icond,:);
    %pre-stim
        gz_pre = [gz_pre, mean(gz_trial(1,1:zeroind))];
    %stim
        gz_stim = [gz_stim, gz_trial(1,zeroind+1:end)];

end
gz_pre = double(gz_pre);
gz_stim = double(gz_stim);
base_pre = double(base_pre);
base_stim = double(base_stim);
%now do the same as did above on averages
figure;
hold on
norespinds = find(base_pre == 0);
p_emerge_pre = size(find(gz_pre(norespinds)~=0),2) ./ size(norespinds,2)
respinds = find(base_pre ~=0);
p_retract_pre = size(find(gz_pre(respinds) == 0),2) / size(respinds,2)
scatter(base_pre(respinds),gz_pre(respinds),30,'b')
scatter(base_pre(norespinds),gz_pre(norespinds),30,'k')
maxrate = (max([max(base_pre) max(gz_pre)]));
cf = fit(base_pre',gz_pre','poly1')
fitline = cf.p1*base_pre+cf.p2;
Rsq = rsquare(gz_pre,fitline);
line(base_pre,fitline,'color','r')
line([0,maxrate],[0,maxrate],'color','k')
set(gca,'YLim',[0,maxrate],'XLim',[0,maxrate])
text(0.5,0.5,['gain = ' num2str(cf.p1) '  R^2 = ' num2str(Rsq)],'FontSize',14)
xlabel('mean pre-stim spike rate baseline','FontSize',14)
ylabel('mean pre-stim spike rate gabazine','FontSize',14)
title('ALL CASES; pre-stimulus period ; trials aligned to stimulus onset','FontSize',14)
legend('pre-stim average','fitline')

%%


 resp_base_stim = base_stim(1,respinds);
resp_gz_stim = gz_stim(1,respinds);
X = [resp_base_stim;resp_gz_stim]';
[IDX,C] = kmeans(X,3);

figure;
hold on
scatter(resp_base_stim(1,find(IDX == 1)),resp_gz_stim(1,find(IDX == 1)),30,'b')
scatter(resp_base_stim(1,find(IDX == 2)),resp_gz_stim(1,find(IDX == 2)),30,'g')
scatter(resp_base_stim(1,find(IDX == 3)),resp_gz_stim(1,find(IDX == 3)),30,'r')


%%

hfig = figure;
hold on
line(Result_base.times,Result_base.FanoFactor,'color','b','LineWidth',3)
line(Result_base_end.times,Result_base_end.FanoFactor,'color','r','LineWidth',3)
legend('Fano base --> onset' ,'Fano stim --> offset')
axis tight
ylims = get(gca,'YLim');
line([0,0],ylims,'color','k')
ylabel('Mean Matched Fano Factor','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

hfig = figure;
hold on
line(Result_base.times,Result_base.FanoFactorAll,'color','b','LineWidth',3)
line(Result_base_end.times,Result_base_end.FanoFactorAll,'color','r','LineWidth',3)
legend('Fano base --> onset' ,'Fano stim --> offset')
axis tight
ylims = get(gca,'YLim');
line([0,0],ylims,'color','k')
ylabel('All Data Fano Factor','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

hfig = figure;
hold on
line(Result_base.times,Result_base.meanRateAll,'color','b','LineWidth',3)
line(Result_base_end.times,Result_base_end.meanRateAll,'color','r','LineWidth',3)
legend('spikerate base --> onset' ,'spikerate stim --> offset')
axis tight
ylims = get(gca,'YLim');
line([0,0],ylims,'color','k')
ylabel('Mean Spike Rate','FontSize',14)
set(hfig,'Position',[388         517    1038         281])

%%
scatterParams.axLim = 'auto';
scatterParams.axLen = 5;
plotScatter(Result_base, -1000, scatterParams); 
text(2.5, 7, '1000 ms before target', 'hori', 'center'); 
plotScatter(Result_base, 500, scatterParams);
text(2.5, 7, '500 ms after target', 'hori', 'center'); 
plotScatter(Result_base, 1500, scatterParams);
text(2.5, 7, '1500 ms after target', 'hori', 'center');



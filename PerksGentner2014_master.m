% Perks and Gentner 2014
% compiled code to generate analyses and figures

% intracellular dataset:  
    % datafolder = 
% extracellular dataset recorded by JVT for Thompson et al 2013: 
    % datafolder = '/Users/kperks/GitHub/iontoncm/'
    %   datafolder contains bird folders each with penetration and site
    %   folders that are called by the script to extract requested spike data
    
r=rigdef('mac');


set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesColorOrder',[0 0 0;1 0 1;0 1 1])
%%
% "repexpts" is a list of experiments for which there were stimuli blocks
% of 5 trials or more in current clamp
repexpts = {
    'KP_B130_131120_p1c1a.mat'
    'KP_B130_131120_p1c2a.mat'
    'KP_B130_131120_p1c2b.mat'
    'KP_B130_131120_p1c2c.mat'
    'KP_B130_131120_p1c3.mat'
    'KP_B136_131205_p1c2.mat'
    'KP_B136_131205_p2c2.mat'
    'KP_B136_131205_p3c1.mat'
    'KP_B580_120824_p1c2_b.mat'
    'KP_B680_131219_p1c1.mat'
    'KP_B680_131219_p1c2.mat'
    'KP_B682_130219_p1c1.mat'
    'KP_B689_131407_p2c1b.mat'
    'KP_B689_131407_p2c1c.mat'
    'KP_B689_131407_p2c1d.mat'
    'KP_B689_131407_p3c1.mat'
    'KP_B694_131212_p1c1.mat'
    'KP_B694_131212_p1c1b.mat'
    'KP_B694_131212_p2c1.mat'
    'KP_B694_131212_p2c2.mat'
    'KP_B694_131212_p3c1.mat'
    'KP_B694_131212_p4c1.mat'
    'KP_B790_140127_p1c1.mat'
    'KP_B790_140127_p1c2.mat'
    'KP_B855_130304_p1c2.mat'};

for icell = 1:size(repexpts,1)
    rootname{icell} = repexpts{icell}(1:19);
    
end
unq_expts = unique(rootname); % some different experiment instances were the same cell, ...
% just different set of stimuli (usually when one set of stimuli were
% longer than the rest because the export script doesn't handle that case
% from IGOR so had to do different experiments for different stimulus
% durations

%% get spikethreshold for each expt and save it in the expt .mat
% RepExptsDat = [];
% for iexpt = 1:size(repexpts,1)
% 
%     load([r.Dir.Expt repexpts{iexpt}])
%     RepExptsDat(iexpt).exptname = expt.name;
%     vmexpt = filtesweeps(expt,0,'Vm',0);
%     hfig = figure;plot(HighpassGeneral(vmexpt.wc.data,1/expt.wc.dt)')
%     RepExptsDat(iexpt).spikethresh = input('spikethresh');
%     close (hfig)
% end
% 
% save('/Users/kperks/GitHub/PerksGentner2014/RepExptsDat_prep.mat','RepExptsDat')

%%
load('/Users/kperks/GitHub/PerksGentner2014/RepExptsDat_prep.mat','RepExptsDat')

%% with repexpts list and unqexpts, make dat struct that will just save and load later
%CreateWholeCellDataset
trials = 5;
clampval = 0;
% make data struct to save in analysis folder for quick use...
% {1,cell}{stimnum,2}
%   .exptname
%   .stimname
%   .dB
%   .data
%   .spikethresh
%   .baselinewin
%   .sigon
%   .siglen
%   .wav
%   .dt

dat = [];

for iunq = 1:size(unqexpts,1)
    
    sigind = 1;
    for icell = 1:size(repexpts,1)
        if strcmp(unqexpts{iunq},repexpts{icell}(1:19))
            load([r.Dir.Expt repexpts{icell}])
            
            %make a list of all stims for which there were more than 5
            %trials
            table=getClampTab(expt,{'clamp',clampval});
            keepsigs=reprequire(table,trials);
            thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
            
            if isempty(thiscond)
                continue
            end
            sigon = round(expt.analysis.params.waveonset_time/expt.wc.dt);
            vmexpt = filtesweeps(expt,0,'Vm',0);
            spikethresh = RepExptsDat(icell).spikethresh;
            
            for istim = 1:size(thiscond,2)
                sigstruct = [];
                sigstruct.exptname = expt.name;
                
%                 strind = regexp(thiscond(istim).wavnames,'d');
%                 sigstruct.stimname = thiscond(thispair(idb)).wavnames(1:strind-1);
                sigstruct.stimname = thiscond(istim).wavnames;
                sigstruct.dB = (thiscond(istim).wavnames(end-1:end));
                sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(istim).wavnames);
                sigstruct.data = sigexpt.wc.data;
                sigstruct.spikethresh = spikethresh;
                siglen = round(size(thiscond(istim).wavs,1)/44100/expt.wc.dt);
                sigstruct.baselinewin = expt.analysis.params.baselinewin;
                sigstruct.sigon = sigon;
                sigstruct.siglen = siglen;
                sigstruct.wav = thiscond(istim).wavs;
                sigstruct.steptime = expt.analysis.params.steptime;
                sigstruct.dt = expt.wc.dt;
                
                
                dat{1,iunq}{sigind} = sigstruct;
                sigind = sigind+1;
            end
            
        end
    end
end

% save('/Users/kperks/GitHub/PerksGentner2014/RepExptsDat_struct.mat','dat')
%%
load('/Users/kperks/GitHub/PerksGentner2014/RepExptsDat_struct.mat','dat')

%% from dat struct get intrinsic and whole cell properties for each cell recorded

PlotRinRsCm = 0;
PlotPop = 0;
plotSpikeShape = 0;
plotResidVm = 0;


%get RinRsCm for each actual unique expt
%would be worth it to go back to the raw data and re-do these
%because raw data has more data point which would help accuracy for Rs
%measurement and therefore Tau/Cm in particular
Rin = [];
Rs = [];
Rv = [];
sagR = [];
TaoCell = [];
TaoV = [];
Cm = [];
expt_stepdata = [];
expt_baselineVm = [];
V_f = [];


for iunq = 1:size(unq_expts,2)
   
    tmpstepdata = [];
    tmpbaselinedata = [];
    this_cell = dat{1,iunq};
    for istim=1:size(this_cell,2)
        
            this_stim = this_cell{istim};
            stepdur = diff(this_stim.steptime); %round(0.25/expt.wc.dt);
            stepstart = 553; %this_stim.steptime(1); %
            
            tmpstepdata = [tmpstepdata;this_stim.data(:,stepstart:stepdur)*1000];
            tmpbaselinedata = [tmpbaselinedata; this_stim.data(:,1:stepstart)*1000];
           
           
        
    end
    
    %only using trials that are in the lower 50% of the std values...
    %%%%%%%%%%%later should make sure that for each experiment there are
    %%%%%%%%%%%enough trials left at that point to get a good estimate...
    %%%%%%%%%%%%%%%%%%could just look at fits to see if reasonable?
    
    %%%%%%%**************
    %%%********DONT SELECT OUT LOWEST STD TRACES BECAUSE UNDERESTIMATING
    %%%HYPERPOL INDUCED CONDUCTANCE**************
    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    %****NO TO THE ABOVE>>> TRIED THAT AND IT IS NOT NOT NOT BETTER
    std_distrib = std(tmpstepdata(:,300:round((stepdur-stepstart)/2))');
    %      figure;hist(std_distrib)
    keepinds = find(std_distrib <= median(std_distrib));
    tmpstepdata = tmpstepdata(keepinds,:);
    expt_stepdata(iunq,:) = mean(tmpstepdata);
    expt_baselineVm(iunq) = mean(mean(tmpbaselinedata));
    
    [out_struct, hfig,hfig2,hfig3] = MetaResponseAnal_RsRin(tmpstepdata,expt);
    allfields = fieldnames(out_struct);
    for ifield = 1:size(allfields,1)
        s = [allfields{ifield} '(iunq) = out_struct.' allfields{ifield} ';'];
        eval(s)
    end
    
    if PlotRinRsCm == 1
        foldername = '/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/RinRsCm/';
        saveas(hfig,[foldername unq_expts{iunq} 'Stepdata.fig'])
        saveas(hfig,[foldername unq_expts{iunq} 'Stepdata.png'])
        
        if ~isempty(hfig2)
            saveas(hfig2,[foldername unq_expts{iunq} 'Stepdata_3sumfit.fig'])
            saveas(hfig2,[foldername unq_expts{iunq} 'Stepdata_3sumfit.png'])
        end
    end
    
    close(hfig)
    close(hfig2)
    close(hfig3)
   
end

for iunq = 1:size(expt_stepdata,1)
    allsteps(iunq,:) = expt_stepdata(iunq,:) - expt_stepdata(iunq,1);
end
figure
line([1:size(allsteps,2)]*expt.wc.dt,allsteps,'color',[0.5 0.5 0.5])

for iunq = 1:size(expt_stepdata,1)
    allsteps(iunq,:) = expt_stepdata(iunq,:) - min(expt_stepdata(iunq,:));
    allsteps(iunq,:) = allsteps(iunq,:) / max(allsteps(iunq,:));
end
figure
line([1:size(allsteps,2)]*expt.wc.dt,allsteps,'color',[0.5 0.5 0.5],'LineWidth',1)
line([1:size(allsteps,2)]*expt.wc.dt,mean(allsteps),'color','k','LineWidth',3)

if PlotPop == 1
    
    [n,p] = empcdf(TaoCell);
    hfig = figure; hold on
    stairs(n*1000,p)
    line([median(TaoCell)*1000,median(TaoCell)*1000],[0,1],'color','r')  %plot in milliseconds (*1000)
    title(['median Tau = ' num2str(median(TaoCell)*1000)])
    xlabel('Passive Membrane Time Constant (ms)')
    ylabel('CDF')
    
    [n,p] = empcdf(TaoV);
    hfig = figure; hold on
    stairs(n*1000,p)
    line([median(TaoV)*1000,median(TaoV)*1000],[0,1],'color','r')  %plot in milliseconds (*1000)
    title(['median Tau = ' num2str(median(TaoV)*1000)])
    xlabel('Voltage-Dependent Membrane Time Constant (ms)')
    ylabel('CDF')
    
    [n,p] = empcdf(Rin);
    hfig = figure; hold on
    stairs(n,p)
    line([median(Rin),median(Rin)],[0,1],'color','r')  %plot in megaOhm (*1000)
    title(['median Rin = ' num2str(median(Rin))])
    xlabel('Input Resistance (megaOhm)')
    ylabel('CDF')
    
    [n,p] = empcdf(Rs);
    hfig = figure; hold on
    stairs(n,p)
    line([median(Rs),median(Rs)],[0,1],'color','r')  %plot in megaOhm (*1000)
    title(['median Rs = ' num2str(median(Rs))])
    xlabel('Series Resistance (megaOhm)')
    ylabel('CDF')
    
    
    figure
    scatter([1:size(unq_expts,2)], Rs)
    set(gca,'XTick',[1:size(unq_expts,2)])
    set(gca,'XTickLabel',unq_expts)
    rotateXLabels(gca,90)
    line([1,20],[60,60],'color','k')
    line([1,20],[50,50],'color','k')
    set(gca,'XLim',[1,20])
    ylabel('access resistance','FontSize',14)
    
    figure;
    scatter(Rin,TaoCell*1000,100,'k','fill')
    ylabel('Time Constant (msec)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)
    
    figure;
    scatter(Rin,Cm,100,'k','fill')
    ylabel('Cell Capacitance (pFarad)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)
    
    figure;
    scatter(Rin,TaoV*1000,100,'k','fill')
    ylabel('Voltage-Dependent Time Constant (msec)','FontSize',14)
    xlabel('Input Resistance (megaOhm)','FontSize',14)  
    
end

%% 
%spike shapes


spk_initVm = []; %this is the output of MetaResponseAnal_Spikes
spk_shape = [];

%restrict this analysis to cells that have a low access resistance
%this will minimize lowpass filtering

for icell=11:size(dat,2)
    if Rs(icell) < 50 % && iexpt ~=22
        this_cell = dat{1,icell};
        for istim = 1:size(this_cell,2)
            this_stim = this_cell{istim};
            sigon = this_stim.sigon;
            sigoff = sigon + this_stim.siglen;
            basetimes = this_stim.baselinewin;
            sigdata = this_stim.data*1000;
            sigdata_filt = medfilt1(sigdata,200,[],2);
            highpassdata=HighpassGeneral(sigdata,1/expt.wc.dt);
            
            Vspkthr = this_stim.spikethresh*1000;
            negative = 0;
            
            %set up input_struct with data to pass analysis functions
            input_struct.dt = this_stim.dt;
            input_struct.spk_thresh = Vspkthr;
            input_struct.sigdata = sigdata(:,basetimes(1):end);
            input_struct.sigdata_filt = sigdata_filt(:,basetimes(1):end);
            input_struct.highpassdata = highpassdata(:,basetimes(1):end);
            
            % get spike shapes and spike threshold and actual spike peak times to get vm during spike
 %%%%%%%% i need to incorporate the new way of doing this into the old analysis
 %%%%%%%%% because right now the spike threshold finding is not as reliable
 %%%%%%%%% as it used to be
            
            [out_struct, hfig] = Anal_Spikes(input_struct);
            if isempty(out_struct)%if there are no spikes this cell:stimulus pair then continue to next
                continue
            end
            
            allfields = fieldnames(out_struct);
            for ifield = 1:size(allfields,1)
                s = [allfields{ifield} '{icell,istim} = out_struct.' allfields{ifield} ';'];
                eval(s)
            end
%             set(hfig,'Visible','on');
            close(hfig)
            
            
        end
    end
end

useSpks = [];
spkwid = [];
spkheight = [];
spkdvdt = [];
%for each cell, collect all spike shapes, etc
for icell = 1:size(spk_shape,1)
    thiscellshape = [];
    thiscellinit = [];
    for istim = 1:size(spk_shape,2)
        thiscellshape = [thiscellshape;spk_shape{icell,istim}];
        thiscellinit = [thiscellinit,spk_initVm{icell,istim}];
    end
    Shape{icell} = thiscellshape;
    Vthresh{icell} = thiscellinit;
    
    
    input_struct = [];
    input_struct.dt = dat{1,icell}{1}.dt;
    input_struct.shape = Shape{icell};
    input_struct.initVm = Vthresh{icell};
    
    [out_struct,hfig] = MetaResponseAnal_SpikeCluster(input_struct);
    
    allfields = fieldnames(out_struct);
    for ifield = 1:size(allfields,1)
        s = [allfields{ifield} '{icell} = out_struct.' allfields{ifield} ';'];
        eval(s)
    end
end

%get mean spike shape and params for each cell averaged over all stimuli unq_spk_shape_norm = [];
unq_spk_shape = [];
unq_wid = [];
unq_height = [];
unq_dvdt = [];
unq_initVm = [];
unq_Vrest = [];
ii = 1;
for icell = 1:size(dat,2)
    
        thisSpks = useSpks(unqinds);
        thiswid = spkwid(unqinds);
        thisheight = spkheight(unqinds);
        thisdvdt = spkdvdt(unqinds);
        thisinitVm = spk_initVm(unqinds);
        thisVrest = VrestSpk(unqinds);
        tmpshape = [];
        tmpwid = [];
        tmpheight = [];
        tmpdvdt = [];
        tmpinitVm = [];
        tmpVrest = [];
        for icond = 1:size(thisSpks,2)
            tmpshape =  [tmpshape; thisSpks{icond}];
            tmpwid =  [tmpwid, thiswid{icond}];
            tmpheight =  [tmpheight, thisheight{icond}];
            tmpdvdt =  [tmpdvdt, thisdvdt{icond}];
            tmpinitVm =  [tmpinitVm, thisinitVm{icond}];
            tmpVrest =  [tmpVrest, thisVrest(icond)];
        end
        
        unq_spk_shape(ii,:) = mean(tmpshape,1);
        unq_spk_shape_norm(ii,:) = mean(tmpshape,1) / max(mean(tmpshape,1));
        unq_wid(ii,:) = mean(tmpwid);
        unq_height(ii,:) = mean(tmpheight);
        unq_dvdt(ii,:) = mean(tmpdvdt);
        unq_initVm(ii,:) = mean(tmpinitVm);
        unq_initVmAll{ii} = tmpinitVm;
        unq_Vrest(ii,:) = mean(tmpVrest);
        ii = ii+1;
    end
end

figure;
line([1:size(unq_spk_shape_norm,2)]*expt.wc.dt*1000,unq_spk_shape_norm,'color',[0.5 0.5 0.5])
line([1:size(unq_spk_shape_norm,2)]*expt.wc.dt*1000,mean(unq_spk_shape_norm),...
    'color','k','LineWidth',5)
axis tight
set(gca,'YLim',[-0.3,1.1])
xlabel('msec','FontSize',14)
ylabel('normalized height','FontSize',14)
title(['n cells = ' num2str(size(unq_spk_shape,1))],'FontSize',14)


figure;
line([1:size(unq_spk_shape,2)]*expt.wc.dt*1000,unq_spk_shape,'color',[0.5 0.5 0.5])
line([1:size(unq_spk_shape,2)]*expt.wc.dt*1000,mean(unq_spk_shape),...
    'color','k','LineWidth',5)
axis tight
%set(gca,'YLim',[-0.3,1.1])
xlabel('msec','FontSize',14)
ylabel('height','FontSize',14)
title(['n cells = ' num2str(size(unq_spk_shape,1))],'FontSize',14)

hfig = figure;
hold on
subplot(3,1,1)
hist(unq_wid)
title('spike width (msec)')
subplot(3,1,2)
hist(unq_height)
title('spike height (mV)')
subplot(3,1,3)
hist(unq_dvdt)
title('dvdt initial rise (V/sec)')
set(hfig,'Position',[98   131   463   663])


X = [unq_wid,unq_height,unq_dvdt];
[idx,ctrs] = kmeans(X,2);

hfig= figure
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),12,'r','fill')
hold on
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),12,'b','fill')
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx',...
    'MarkerSize',12,'LineWidth',2)
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'ko',...
    'MarkerSize',12,'LineWidth',2)
xlabel('spkwid msec')
ylabel('spkheight')
zlabel('spkdvdt V/sec')
title(['n = ' num2str(size(idx,1))]);

%% response window analysis on membrane potential and spiking

% (can load this:
% load([r.Dir.Expt 'Analysis/Meta_Responses_DB/exptsForICanal_revised2.mat'])
% to look at this analysis without re-doing the prep for it)

% right now, this analysis uses clipping spikes method instead of median
% filter to get rid of spikes

% i had picked out slightly different spike thresholds before so i have
% that information here in case i need to go back an change the ".dat"
% dataset to use the spikethresholds used in the paper
% these match up to "repexpts"
spkthreshAll = [0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.01
0.025
0.01
0.01
0.01
0.01
0.01
0.01
0.015
0.01
0.01
0.01
0.01
0.005
0.01]*1000;
%%
doforEPS = 0;
doplot = 0;
Vthresh = [];
BaseVmVar = [];
VmRespVec_up = [];
VmRespVec_low = [];
VmBaseConf = [];
SubRespWin_up = [];
SubRespWin_low = [];
SpkRespWin_up = [];
SpkRespWin_low = [];
SpkRespRate_up = [];
SpkRespRate_low = [];
SpkRespRate_not = [];
VmRespMean_up = [];
VmRespVar_up = [];
VmRespMean_low = [];
VmRespVar_low = [];
VmNotVec = [];
VmNotVar = [];
siglen = [];
signame = [];

TrackCellInfo = [];
TrackSigInfo = [];
%%%%need to keep track of stimulus names so that at end can only report for
%%%%each stimulus once (pick the loudest db if there are multiple)
%% this analysis still uses "repexpts" and "unqexpts"
%%%%%%%%%%% switch it over to pulling from dat.mat dataset struct
unqind = 1;
for iunq = 1:size(unq_expts,2)
    unq_expts{iunq}
    
    allstepdata = [];
    this_cell = unq_expts{iunq};
    Vrest = [];
    response_vm_confint = [];
    response_vm_mean = [];
    response_vm_max = [];
    up_win = [];
    low_win = [];
    allVmRespVec_up = [];
    allVmRespVec_low = [];
    
    sigind = 1;
    for iexpt=1:size(repexpts,1)
        if ~isempty(regexp(repexpts{iexpt},this_cell))
            thisexpt=repexpts{iexpt};
            load([r.Dir.Expt thisexpt])
                        
            vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
            % for now do not care if db or not. else: [dbstimcond,dblevels]=getDBstimcond(vmexpt);
            table=getClampTab(expt,{'clamp',0});
            keepsigs=reprequire(table,trials);
            allsig=table.sigsplayed;
            repsigs=allsig(keepsigs);
            stimcond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
            
            %get spiking threshold
            spk_thresh = spkthreshAll(iexpt);
            
            for istim = 1:size(stimcond,2)
                input_struct = [];
                %skip warped shortened tempo stims
                if ~isempty(regexp(stimcond(istim).wavnames,'ws'))
                    continue
                end
                 TrackCellInfo{iunq,sigind} = expt.name;
                  TrackSigInfo{iunq,sigind} = stimcond(istim).wavnames;
                  
                sigexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(istim).wavnames);
                sigdata = sigexpt.wc.data(:,1:end-4)*1000;
                sigdata_filt = medfilt1(sigdata,200,[],2);
                
                input_struct.sigdata = sigdata;
                input_struct.sigdata_filt = sigdata_filt;
                
                
                basetimes = sigexpt.analysis.params.baselinewin;
                [sigon,sigoff]=GetSigTimes(sigexpt,stimcond,istim);
                
                %set up input_struct with data to pass analysis functions
                input_struct.basetimes = basetimes;
                input_struct.sigon = sigon;
                input_struct.sigoff = sigoff;
                
                siglen{unqind,sigind} = sigoff-sigon; %size(sigdata_filt(:,basetimes(1):end),2);
                signame{unqind,sigind} = stimcond(istim).wavnames;
                
                %prep spiking to be continuous and get spiking windows
                highpassdata=HighpassGeneral(sigdata,1/expt.wc.dt);
                [spikesmat, gausstosmooth, spiketimes]=getspikesmat(highpassdata,spk_thresh,expt.wc.dt);
                
                if ~isempty(gausstosmooth)
                    spkvec = zeros(size(spikesmat,1),size(spikesmat,2));
                    spkvec(find(spikesmat)) = 1;
                    smoothspk = [];
                    for itrial = 1:size(spkvec,1)
                        smoothspk(itrial,:) = conv(spkvec(itrial,:),gausstosmooth,'same');
                    end
                    
                    conf_p = 95;
                    windowsize = 500;
                    binsize = 10;
                    bin_p = 0.85;
                    
                    confint_spk = getCDFconf (mean(smoothspk(:,basetimes(1):sigon)),conf_p);
                    
                    [up_inds_spk, low_inds_spk] = WindowResponse(mean(smoothspk(:,sigon:sigoff)),...
                        confint_spk, windowsize, binsize, bin_p);
                    all_inds_spk = union(up_inds_spk,low_inds_spk);
                    all_inds = [1:(sigoff-sigon)];
                    all_inds(all_inds_spk) = 0;
                    notinds = find(all_inds);
                    up_win_spk = getWindowEdges (up_inds_spk, 1, 1)+sigon;
                    low_win_spk = getWindowEdges (low_inds_spk, 1, 1)+sigon;
                    SpkRespWin_up{unqind,sigind} = up_win_spk;
                    SpkRespWin_low{unqind,sigind} = low_win_spk;
                    
                    %get spike rate per response
                    allrate = [];
                    for iwin = 1:size(up_win_spk,2)
                        allspk = size(find(spikesmat(:,up_win_spk(1,iwin):up_win_spk(2,iwin))==1),1);
                        allrate(iwin) = (allspk/size(spikesmat,1)) / ...
                            (diff(up_win_spk(:,iwin))*expt.wc.dt);
                    end
                    SpkRespRate_up{unqind,sigind} = allrate;
                    
                    allrate = [];
                    for iwin = 1:size(low_win_spk,2)
                        allspk = size(find(spikesmat(:,low_win_spk(1,iwin):low_win_spk(2,iwin))==1),1);
                        allrate(iwin) = (allspk/size(spikesmat,1)) / ...
                            (diff(low_win_spk(:,iwin))*expt.wc.dt);
                    end
                    SpkRespRate_low{unqind,sigind} = allrate;
                    
                    tmp_mat = spikesmat;
                    tmp_mat(:,all_inds_spk) = 0;
                    %double check that sigon:sigoff was not already
                    %accounted for
                    tmp_mat = tmp_mat(:,sigon:sigoff);
                    not_spks = size(find(tmp_mat == 1),1);
                    SpkRespRate_not{unqind,sigind} = (not_spks/size(spikesmat,1))/((sigoff-sigon)*expt.wc.dt);
                    
                    
                    
                    
                    if doplot == 1
                        spkrespfig = figure;
                        hold on
                        for itrial = 1:size(smoothspk,1)
                            thesespks = find(spkvec(itrial,:))*expt.wc.dt
                            for ispike = 1:size(thesespks,2)
                                line([thesespks(ispike),thesespks(ispike)],...
                                    [1*itrial,(1*itrial)+1],'color','k')
                            end
                        end
                        set(gca,'YLim',[1,size(smoothspk,1)+1],'XLim',[1*expt.wc.dt,size(smoothspk,2)*expt.wc.dt])
                        
                        SigTimeBox(gca, (sigon)*expt.wc.dt,sigoff*expt.wc.dt, get(gca,'YLim'),[0.5 0.5 0.5]);
                        for iresp=1:size(up_win_spk,2)
                            SigTimeBox(gca, up_win_spk(1,iresp)*expt.wc.dt, ...
                                up_win_spk(2,iresp)*expt.wc.dt, get(gca,'YLim'),'r');
                        end
                        for inhib=1:size(low_win_spk,2)
                            SigTimeBox(gca, low_win_spk(1,inhib)*expt.wc.dt, ...
                                low_win_spk(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
                        end
                        %                     axis tight
                        set(gca,'TickDir','out')
                        
                        box off
                        set(spkrespfig,'Position',[212         523        1168         283])
                        title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                        
                        %                 saveas(spkrespfig, ...
                        %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/SpkRespWin/' ...
                        %                     expt.name '_' stimcond(istim).wavnames '.fig']);
                        %                 saveas(spkrespfig, ...
                        %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/SpkRespWin/' ...
                        %                     expt.name '_' stimcond(istim).wavnames '.tif']);
                        close(spkrespfig)
                        %
                        %plot without patches for eps
                        if doforEPS == 1;
                            spkrespfig = figure;
                            hold on
                            for itrial = 1:size(smoothspk,1)
                                thesespks = find(spkvec(itrial,:))*expt.wc.dt
                                for ispike = 1:size(thesespks,2)
                                    line([thesespks(ispike),thesespks(ispike)],[1*itrial,(1*itrial)+1])
                                end
                            end
                            set(gca, 'YLim',[-0.5,size(smoothspk,1)+1],'XLim',[1*expt.wc.dt,size(smoothspk,2)*expt.wc.dt])
                            
                            line([(sigon)*expt.wc.dt,sigoff*expt.wc.dt], [1,1],'color','k');
                            for iresp=1:size(up_win_spk,2)
                                line([up_win_spk(1,iresp)*expt.wc.dt, ...
                                    up_win_spk(2,iresp)*expt.wc.dt], [0.5,0.5],'color','r');
                            end
                            for inhib=1:size(low_win_spk,2)
                                line([low_win_spk(1,inhib)*expt.wc.dt, ...
                                    low_win_spk(2,inhib)*expt.wc.dt], [0,0],'color','b');
                            end
                            box off
                            set(gca,'TickDir','out')
                            set(spkrespfig,'Position',[212         523        1168         283])
                            title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                            
                            %                         saveas(spkrespfig, ...
                            %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/SpkRespWin/foreps' ...
                            %                     expt.name '_' stimcond(istim).wavnames '.fig']);
                            %                 saveas(spkrespfig, ...
                            %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/SpkRespWin/foreps' ...
                            %                     expt.name '_' stimcond(istim).wavnames '.tif']);
                            close(spkrespfig)
                        end
                    end
                end
                if isempty(gausstosmooth) %if there were not enough spikes, then there are no spike response windows
                    up_win_spk = [];
                    low_win_spk = [];
                    SpkRespWin_up{unqind,sigind} = up_win_spk;
                    SpkRespWin_low{unqind,sigind} = low_win_spk;
                end
                
                %
                cutdata = [];
                vthresh_sig = [];
                spkind = 1;
                for itrial=1:size(spikesmat,1)
                    spks_trial = spiketimes{itrial};
                    thistrial = sigdata(itrial,:);
                    for ispike = 1:size(spks_trial,2)
                        t1 = spks_trial(ispike);
                        
                        spkwin_size = round((20/1000/dt)/2);
                        if ispike == 1
                            altbegin = [(t1 - spkwin_size),1];
                        else altbegin = [(t1 - spkwin_size),1,spks_trial(ispike-1)+10];
                        end
                        
                        if ispike == size(spks_trial,2);
                            altend = [(t1 + spkwin_size),size(highpassdata,2),];
                        else altend = [(t1 + spkwin_size),size(highpassdata,2),spks_trial(ispike+1)];
                        end
                        
                        spk_win = [max(altbegin),min(altend)];
                        
                        tmpshape = thistrial(spk_win(1):spk_win(2));
                        spkt = t1-spk_win(1);
                        peakt = min(find(tmpshape==max(tmpshape(spkt-10:spkt+20))));
                        dvdt = diff(tmpshape(1:peakt));
                        dvdt_max = max(dvdt(1,end-10:end));
                        peakind = max(find(dvdt == dvdt_max));
                        %         dvdt_thr = 0.033 * dvdt_max;
                        dvdt_thr = 0.1 * dvdt_max;
                        
                        spk_init = max(find(dvdt(1:peakind)<=dvdt_thr));
                        
                        if isempty(spk_init)
                            continue
                        end
                        
                        vthresh_sig(spkind) = tmpshape(spk_init);
                        spk_end = min(find(tmpshape(peakt:end)<...
                            tmpshape(spk_init)))+peakt-1;
                        %if no value less than spk_init in that window)
                        if isempty(spk_end)
                            spk_end = min(find(tmpshape(peakt:end) == ...
                                min(tmpshape(peakt:end)))) + peakt-1;
                        end
                        
                        x = [spk_init,spk_end];
                        xi = [spk_init:1:spk_end];
                        y = [tmpshape(spk_init),tmpshape(spk_end)];
                        yi = interp1(x,y,xi);
                        
                        cutshape = [tmpshape(1:spk_init-1),yi,tmpshape((spk_end+1):end)];
                        
                        startcut = spk_init + spk_win(1) -2;
                        stopcut = spk_end + spk_win(1);
                        thistrial = [thistrial(1:startcut),yi,sigdata(itrial,stopcut:end)];
                        
                        spkind = spkind + 1;
                    end
                    cutdata(itrial,:) = thistrial;
                end
                    Vthresh{unqind,sigind} = vthresh_sig;
                    
                    
                    % get vm response windows
                    
                    out_struct = MetaResponseAnal_VmResponseWin(expt,input_struct);
                    allfields = fieldnames(out_struct);
                    for ifield = 1:size(allfields,1)
                        s = [allfields{ifield} '{sigind} = out_struct.' allfields{ifield} ';'];
                        eval(s)
                    end
                    
                    Vrest{sigind} = mean(min(input_struct.sigdata_filt(:,basetimes(1):basetimes(2))'));
                    
                    VmBaseConf{unqind,sigind} = response_vm_confint{sigind};
                    BaseVmVar{unqind,sigind} = spont_var{sigind};
                    VrestCell{unqind,sigind} = Vrest{sigind};
                    SubRespWin_up{unqind,sigind} = up_win{sigind};
                    SubRespWin_low{unqind,sigind} = low_win{sigind};
                    VmRespVec_up{unqind,sigind} = allVmRespVec_up{sigind};
                    VmRespVec_low{unqind,sigind} = allVmRespVec_low{sigind};
                    VmRespMean_up{unqind,sigind} = vm_mean_up{sigind};
                    VmRespVar_up{unqind,sigind} = vm_var_up{sigind};
                    VmRespMean_low{unqind,sigind} = vm_mean_low{sigind};
                    VmRespVar_low{unqind,sigind} = vm_var_low{sigind};
                    VmNotVec{unqind,sigind} = allVmNotVec{sigind};
                    VmNotVar{unqind,sigind} = allVmNotVar{sigind};
                    
                    if doplot ==1
                        
                        respfig = figure;
                        hold on
                        scaleticks = 1;
                        ydatabound = [min(mean(input_struct.sigdata_filt(:,basetimes(1):end))), max(mean(input_struct.sigdata_filt(:,basetimes(1):end)))];
                        xtime=[1:size(input_struct.sigdata_filt,2)]*expt.wc.dt;
                        line(xtime,mean(input_struct.sigdata_filt),'color','k','LineWidth',3);
                        plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(2),response_vm_confint{sigind}(2)],'--','color','k')
                        plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(1),response_vm_confint{sigind}(1)],'--','color','k')
                        SigTimeBox(gca, (sigon)*expt.wc.dt,sigoff*expt.wc.dt, get(gca,'YLim'),[0.5 0.5 0.5]);
                        for iresp=1:size(up_win{sigind},2)
                            SigTimeBox(gca, up_win{sigind}(1,iresp)*expt.wc.dt, ...
                                up_win{sigind}(2,iresp)*expt.wc.dt, get(gca,'YLim'),'r');
                        end
                        for inhib=1:size(low_win{sigind},2)
                            SigTimeBox(gca, low_win{sigind}(1,inhib)*expt.wc.dt, ...
                                low_win{sigind}(2,inhib)*expt.wc.dt, get(gca,'YLim'),'b');
                        end
                        axis tight
                        ylims = [Vrest{sigind},ydatabound(2)];
                        xlims = get(gca,'XLim');
                        xlims = [basetimes(1)*expt.wc.dt,xlims(2)];
                        set(gca,'XLim',xlims);
                        
                        set(gca,'YLim',ylims,'YTick',...
                            [(floor(Vrest{sigind})-mod(floor(Vrest{sigind}),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
                            'XTick',[0:1:floor(xtime(end))],'TickDir','out')
                        text(xtime(1,basetimes(1)),round(Vrest{sigind}),[num2str(round(Vrest{sigind})) 'mV'],...
                            'HorizontalAlignment','center',	'BackgroundColor', 'k',...
                            'color',[1,1,1]);
                        box off
                        set(respfig,'Position',[212         523        1168         283])
                        title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                        
                        %                 saveas(respfig, ...
                        %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/' ...
                        %                     expt.name '_' stimcond(istim).wavnames '.fig']);
                        %                 saveas(respfig, ...
                        %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/' ...
                        %                     expt.name '_' stimcond(istim).wavnames '.tif']);
                        close(respfig)
                        if doforEPS == 1;
                            respfig = figure;
                            hold on
                            scaleticks = 1;
                            ydatabound = [min(mean(input_struct.sigdata_filt(:,basetimes(1):end))), max(mean(input_struct.sigdata_filt(:,basetimes(1):end)))];
                            xtime=[1:size(input_struct.sigdata_filt,2)]*expt.wc.dt;
                            line(xtime,mean(input_struct.sigdata_filt),'color','k','LineWidth',3);
                            plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(2),response_vm_confint{sigind}(2)],'--','color','k')
                            plot([xtime(1),xtime(end)],[response_vm_confint{sigind}(1),response_vm_confint{sigind}(1)],'--','color','k')
                            line([ (sigon)*expt.wc.dt,sigoff*expt.wc.dt], [ydatabound(1)-1, ydatabound(1)-1],...
                                'color','k','LineWidth',3);
                            for iresp=1:size(up_win{sigind},2)
                                line([up_win{sigind}(1,iresp)*expt.wc.dt, ...
                                    up_win{sigind}(2,iresp)*expt.wc.dt],[ydatabound(1)-1.2, ydatabound(1)-1.2],...
                                    'color' ,'r','LineWidth',0.5);
                            end
                            for inhib=1:size(low_win{sigind},2)
                                line([low_win{sigind}(1,inhib)*expt.wc.dt, ...
                                    low_win{sigind}(2,inhib)*expt.wc.dt],[ydatabound(1)-1.4, ydatabound(1)-1.4],...
                                    'color' ,'b','LineWidth',0.5);
                            end
                            
                            axis tight
                            Vrest{sigind} = mean(min(input_struct.sigdata_filt(:,basetimes(1):basetimes(2))'));
                            ylims = [ydatabound(1)-2,ydatabound(2)];
                            xlims = get(gca,'XLim');
                            xlims = [basetimes(1)*expt.wc.dt,xlims(2)];
                            set(gca,'XLim',xlims);
                            
                            set(gca,'YLim',ylims,'YTick',...
                                [(floor(Vrest{sigind})-mod(floor(Vrest{sigind}),5)):5:(ceil(ymax)+mod(ceil(ymax),5))],...
                                'XTick',[0:1:floor(xtime(end))],'TickDir','out')
                            text(xtime(1,basetimes(1)),round(Vrest{sigind}),[num2str(round(Vrest{sigind})) 'mV'],...
                                'HorizontalAlignment','center',	'BackgroundColor', 'k',...
                                'color',[1,1,1]);
                            box off
                            set(respfig,'Position',[212         523        1168         283])
                            
                            title([expt.name stimcond(istim).wavnames],'Interpreter','none')
                            %
                            %                 saveas(respfig, ...
                            %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/foreps_' ...
                            %                     expt.name '_' stimcond(istim).wavnames '.fig']);
                            %                 saveas(respfig, ...
                            %                     ['/Users/kperks/GitHub/Data_Mat/Analysis/Meta_Responses_DB/VmRespWin/foreps_' ...
                            %                     expt.name '_' stimcond(istim).wavnames '.tif']);
                            close(respfig)
                        end
                    end
                    sigind = sigind + 1;
                end
            end
        end
        
        %summary data for each expt (cell/signal block)
        %VmVector of all Vm values during responses (up or down)
        
        
        
        unqind = unqind+1;
        
    
end
%% get statts from Vm and spiking response window structs on response duration and magnitude
%%%%%%% and plot them
% Vthresh
% SpkRespWin_up
% SpkRespWin_low
% SpkRespRate_up
% SpkRespRate_low
% SpkRespRate_not
% VmBaseConf
% VrestCell
% SubRespWin_up
% SubRespWin_low
% VmRespVec_up
% VmRespVec_low
% siglen
% signame
% VmRespMean_up
% VmRespVar_up
% VmRespMean_low
% VmRespVar_low
% VmNotVec
% VmNotVar
% BaseVmVar

%distribution of response epoch durations
%mean across all epochs across all cells, not mean of mean per cell

%%%%%%%%%%%%%%% response increases
vm_alluplen = [];
sp_alluplen = [];
flagSpkResp = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_alluplen = [vm_alluplen, diff(SubRespWin_up{iunq,isig})];
            sp_alluplen = [sp_alluplen, diff(SpkRespWin_up{iunq,isig})];
            if ~isempty(find(diff(SpkRespWin_up{iunq,isig})*expt.wc.dt>3))
                flagSpkResp(iunq,isig) = 1;
            end
        end
    end
end
sp_meduplen = median(sp_alluplen) * expt.wc.dt
sp_ci = getCDFconf(sp_alluplen*expt.wc.dt,95)
vm_meduplen = median(vm_alluplen) * expt.wc.dt
vm_ci = getCDFconf(vm_alluplen*expt.wc.dt,95)
figure; hold on
[x,p] = empcdf(vm_alluplen* expt.wc.dt);
stairs(x,p,'color','r')
[x,p] = empcdf(sp_alluplen* expt.wc.dt);
stairs(x,p,'color','k')
legend({'vm "up" epoch durations', 'spk "up" epoch durations'})
xlabel('seconds')

%percent of the stimulus that is facilatory responses
vm_percent_up = [];
sp_percent_up = [];
flagSpkResp = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_percent_up = [vm_percent_up, sum(diff(SubRespWin_up{iunq,isig}))/siglen{iunq,isig}];
            sp_percent_up = [sp_percent_up, sum(diff(SpkRespWin_up{iunq,isig}))/siglen{iunq,isig}];
        end
    end
end
sp_medpercup = median(sp_percent_up) * 100
sp_ci = getCDFconf(sp_percent_up*100,95)
vm_medpercup = median(vm_percent_up) * 100
vm_ci = getCDFconf(vm_percent_up*100,95)
[h,p] = kstest2(vm_percent_up,sp_percent_up)
numberstim = size(sp_percent_up)
numberstim = size(vm_percent_up)

%percent of the stimulus that is suppressive responses
vm_percent_dn = [];
sp_percent_dn = [];
flagSpkResp = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_percent_dn = [vm_percent_dn, sum(diff(SubRespWin_low{iunq,isig}))/siglen{iunq,isig}];
            sp_percent_dn = [sp_percent_dn, sum(diff(SpkRespWin_low{iunq,isig}))/siglen{iunq,isig}];
        end
    end
end
sp_medpercdn = median(sp_percent_dn) * 100
sp_ci = getCDFconf(sp_percent_dn*100,95)
vm_medpercdn = median(vm_percent_dn) * 100
vm_ci = getCDFconf(vm_percent_dn*100,95)
[h,p] = kstest2(vm_percent_dn,sp_percent_dn)
numberstim = size(sp_percent_dn)
numberstim = size(vm_percent_dn)

%percent of the stimulus that is not driving subthreshold responses
vm_percent_not = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_percent_not = [vm_percent_not, size(VmNotVec{iunq,isig},2)/siglen{iunq,isig}];
        end
    end
end
vm_medpercnot = median(vm_percent_not) * 100
vm_ci = getCDFconf(vm_percent_not*100,95)
numberstim = size(vm_percent_not)



% distribution of facilatory response magnitudes
%mean across all epochs across all cells, not mean of mean per cell
vm_allupmag = [];
sp_allupmag = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_allupmag = [vm_allupmag, (VmRespMean_up{iunq,isig})-VrestCell{iunq,isig}];
            sp_allupmag = [sp_allupmag, (SpkRespRate_up{iunq,isig})];
        end
    end
end
sp_medupmag = median(sp_allupmag) 
sp_ci = getCDFconf(sp_allupmag,95)
vm_medupmag = median(vm_allupmag) 
vm_ci = getCDFconf(vm_allupmag,95)
figure; hold on
[x,p] = empcdf(vm_allupmag); 
stairs(x,p,'color','r')
[x,p] = empcdf(sp_allupmag);
stairs(x,p,'color','k')
legend({'vm "up" mean resp mV', 'spk "up" mean resp spk/sec'})

%%%%%%%%%%% response decreases
vm_alldnlen = [];
sp_alldnlen = [];
flagSpkResp = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_alldnlen = [vm_alldnlen, diff(SubRespWin_low{iunq,isig})];
            sp_alldnlen = [sp_alldnlen, diff(SpkRespWin_low{iunq,isig})];
            if ~isempty(find(diff(SpkRespWin_up{iunq,isig})*expt.wc.dt>3))
                flagSpkResp(iunq,isig) = 1;
            end
        end
    end
end
sp_meddnlen = median(sp_alldnlen) * expt.wc.dt
sp_ci = getCDFconf(sp_alldnlen*expt.wc.dt,95)
vm_meddnlen = median(vm_alldnlen) * expt.wc.dt
vm_ci = getCDFconf(vm_alldnlen*expt.wc.dt,95)
figure; hold on
[x,p] = empcdf(vm_alldnlen* expt.wc.dt);
stairs(x,p,'color','r')
[x,p] = empcdf(sp_alldnlen* expt.wc.dt);
stairs(x,p,'color','k')
legend({'vm "low" epoch durations', 'spk "low" epoch durations'})
xlabel('seconds')

% distribution of response magnitudes
%mean across all epochs across all cells, not mean of mean per cell
vm_alldnmag = [];
sp_alldnmag = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            vm_alldnmag = [vm_alldnmag, (VmRespMean_low{iunq,isig})-VrestCell{iunq,isig}];
            sp_alldnmag = [sp_alldnmag, (SpkRespRate_low{iunq,isig})];
        end
    end
end
sp_meddnmag = median(sp_alldnmag) 
sp_ci = getCDFconf(sp_alldnmag,95)
vm_meddnmag = median(vm_alldnmag) 
vm_ci = getCDFconf(vm_alldnmag,95)
figure; hold on
[x,p] = empcdf(vm_alldnmag); 
stairs(x,p,'color','r')
[x,p] = empcdf(sp_alldnmag);
stairs(x,p,'color','k')
legend({'vm "low" mean resp mV', 'spk "low" mean resp spk/sec'})

allspkthresh = [];
for iunq = 1:size(signame,1)
    x_depol = [];
    y_var = [];
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig})
            if ~isempty(Vthresh{iunq,isig})
                allspkthresh = [allspkthresh,-VrestCell{iunq,isig} + Vthresh{iunq,isig}];
            end
        end
    end
end
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%% group of anal for covariance of variance and repsonse
%%%%%%%%%%%% magnitude


% vm_allupmag = [];
% vm_allupdepol = [];
% vm_allupvar = [];
figure;hold on

rsq = [];
rmse = [];
scale = [];
all_depol = [];

% sp_allupmag = [];  %%%%%need to do calc of spiking responses in previous script above
for iunq = 1:size(signame,1)
    x_depol = [];
    y_var = [];
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});
            if ~isempty(VmRespMean_up{iunq,isig})
                %         vm_allupmag = [vm_allupmag, VmRespMean_up{iunq,isig}];
                %         vm_allupdepol = [vm_allupdepol, repmat(-VrestCell{iunq,isig},1,size(VmRespMean_up{iunq,isig},2)) ...
                %             + VmRespMean_up{iunq,isig}];
                %         vm_allupvar = [vm_allupvar, VmRespVar_up{iunq,isig}];
                x_depol = [x_depol,(repmat(-VrestCell{iunq,isig},1,size(VmRespMean_up{iunq,isig},2)) ...
                    + VmRespMean_up{iunq,isig})];
                all_depol = [all_depol, (repmat(-VrestCell{iunq,isig},1,size(VmRespMean_up{iunq,isig},2)) ...
                    + VmRespMean_up{iunq,isig})];
                y_var = [y_var,(VmRespVar_up{iunq,isig}./BaseVmVar{iunq,isig})];
                
            end
        end
    end
    up_vm_var(iunq) = mean(y_var);
    
    scatter(x_depol,y_var)
    
    cf = fit(x_depol',y_var','poly1');
    % fitline = cf.p1*dnsampRMS + cf.p2;
    fitline = cf.p1*x_depol + cf.p2;
    line(x_depol,fitline,'color','r')
    scale(iunq) = cf.p1;
    [rsq(iunq) rmse(iunq)] = rsquare(y_var,fitline);    
end
set(gca,'FontSize',16)
title(' "depolarizing" responses')
ylabel('proportion spontaneous variance')
xlabel('mean response magnitude (mV above Vrest)')
% vm_medupmag = median(vm_allupmag)
% vm_ci = getCDFconf(vm_allupmag,95)
% figure; hold on
% [x,p] = empcdf(vm_allupmag);
% stairs(x,p,'color','r')
% vm_medupdepol = median(vm_allupdepol)
% vm_ci = getCDFconf(vm_allupdepol,95)
% figure;
% scatter(vm_allupdepol,vm_allupvar)
figure;hold on
   
    rsq = nan(size(signame,1),size(signame,2));
    rmse = nan(size(signame,1),size(signame,2));
    scale = nan(size(signame,1),size(signame,2));
    all_hyper = [];
% sp_allupmag = [];  %%%%%need to do calc of spiking responses in previous script above
for iunq = 1:size(signame,1)
 x_depol = [];
    y_var = [];
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});
            if ~isempty(VmRespMean_low{iunq,isig})
                if (VmRespMean_low{iunq,isig}) ~= 0
                %         vm_allupmag = [vm_allupmag, VmRespMean_up{iunq,isig}];
                %         vm_allupdepol = [vm_allupdepol, repmat(-VrestCell{iunq,isig},1,size(VmRespMean_up{iunq,isig},2)) ...
                %             + VmRespMean_up{iunq,isig}];
                %         vm_allupvar = [vm_allupvar, VmRespVar_up{iunq,isig}];
                x_depol = [x_depol,(repmat(-VrestCell{iunq,isig},1,size(VmRespMean_low{iunq,isig},2)) ...
                    + VmRespMean_low{iunq,isig})];
                all_hyper = [all_hyper,(repmat(-VrestCell{iunq,isig},1,size(VmRespMean_low{iunq,isig},2)) ...
                    + VmRespMean_low{iunq,isig})];
                y_var = [y_var,(VmRespVar_low{iunq,isig}./BaseVmVar{iunq,isig})];
                end
            end
        end
    end
    low_vm_var(iunq) = mean(y_var);
    
    if size(x_depol,2) > 1
    scatter(x_depol,y_var)

        cf = fit(x_depol',y_var','poly1');
        % fitline = cf.p1*dnsampRMS + cf.p2;
        fitline = cf.p1*x_depol + cf.p2;
        line(x_depol,fitline,'color','r')
        scale(iunq) = cf.p1;
        [rsq(iunq) rmse(iunq)] = rsquare(y_var,fitline);
    end
   
end
set(gca,'FontSize',16)
title(' "hyperpolarizing" responses')
ylabel('proportion spontaneous variance')
xlabel('mean response magnitude (mV above Vrest)')

figure;hold on
rsq = nan(1,size(signame,1));
rmse = nan(1,size(signame,1));
scale = nan(1,size(signame,1));
all_not = [];
% sp_allupmag = [];  %%%%%need to do calc of spiking responses in previous script above
for iunq = 1:size(signame,1)
    x_depol = [];
    y_var = [];
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});
            if ~isempty(VmNotVar{iunq,isig})
                %         vm_allupmag = [vm_allupmag, VmRespMean_up{iunq,isig}];
                %         vm_allupdepol = [vm_allupdepol, repmat(-VrestCell{iunq,isig},1,size(VmRespMean_up{iunq,isig},2)) ...
                %             + VmRespMean_up{iunq,isig}];
                %         vm_allupvar = [vm_allupvar, VmRespVar_up{iunq,isig}];
                x_depol = [x_depol,(-VrestCell{iunq,isig} ...
                    + mean(VmNotVec{iunq,isig}))];
                all_not = [all_not,(-VrestCell{iunq,isig} ...
                    + mean(VmNotVec{iunq,isig}))];
                y_var = [y_var,(mean(VmNotVar{iunq,isig})./BaseVmVar{iunq,isig})];
                
            end
        end
    end
    not_vm_var(iunq) = mean(y_var);
    
    if size(x_depol,2) > 1
        scatter(x_depol,y_var)
        
        cf = fit(x_depol',y_var','poly1');
        % fitline = cf.p1*dnsampRMS + cf.p2;
        fitline = cf.p1*x_depol + cf.p2;
        line(x_depol,fitline,'color','r')
        scale(iunq) = cf.p1;
        [rsq(iunq) rmse(iunq)] = rsquare(y_var,fitline);
    end
end
set(gca,'FontSize',16)
title(' "not" responses')
ylabel('proportion spontaneous variance')
xlabel('mean response magnitude (mV above Vrest)')


vm_notvar = [];
vm_notdepol = [];
for iunq = 1:size(signame,1)
    for isig = 1:size(signame,2) %the array is size of max number stimuli... need to go through and ask if empty for each
        if ~isempty(signame{iunq,isig});   
            if ~isempty(VmNotVec{iunq,isig})
        vm_notdepol = [vm_notdepol, -VrestCell{iunq,isig}...
            + median(VmNotVec{iunq,isig})];
        vm_notvar = [vm_notvar, median(VmNotVar{iunq,isig})];
            end
        end
    end
end
vm_mednotdepol = median(vm_notdepol)
vm_ci = getCDFconf(vm_notdepol,95)
figure;
scatter(vm_notdepol,vm_notvar)
ylabel('var of not repsonses')
xlabel('vm of not responses')


figure; hold on
scatter(repmat(1,1,size(low_vm_var,2)),low_vm_var)
scatter(repmat(2,1,size(not_vm_var,2)),not_vm_var)
scatter(repmat(3,1,size(up_vm_var,2)),up_vm_var)
plot([0,4],[1,1],'--k')
set(gca,'XTick',[1:3],'XTickLabel',{'hyper','spont','depol'},...
    'XLim',[0.5,3.5],'FontSize',16)
ylabel('proportion spontaneous variance')

figure; hold on
edges = [-50:2:50];
x_depol = histc(all_depol,edges)/size(all_depol,2)*100;
x_hyper = histc(all_hyper,edges)/size(all_hyper,2)*100;
x_not = histc(all_not,edges)/size(all_not,2)*100;
x_spk = histc(allspkthresh,edges)/size(allspkthresh,2)*100;
spkbins = find(x_spk>0.1);
bar(edges(1,spkbins)+1,x_depol(1,spkbins),'y','BarWidth',1)
stairs(edges,x_depol,'color','r','LineWidth',3)
stairs(edges,x_hyper,'color','b','LineWidth',3)
stairs(edges,x_not,'color','k','LineWidth',3)
 stairs(edges,x_spk,'color',[0.5 0.5 0.5],'LineWidth',3)
set(gca,'XLim',[-5,40],'YLim',[0,100],'FontSize',16)
legend({[num2str(sum(x_depol(1,spkbins))) '% overlap'],'depol','hyper','not','Vthresh'})
xlabel('mean magnitude of response (mV above Vrest)')
ylabel('proportion total responses')

depolVmOverlapSpk = sum(x_depol(1,spkbins))
CI_97_spkTHresh = getCDFconf(allspkthresh,97)


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%% end of this group of anal for covariance of variance and
%%%%%%%%%%% response


  
%% K_means clustering test of membrane voltage trace specificity
% original .m file is: 'KmeansUniqueness.m'
% data from this saves as :
% matnames = {'KmeansAccSeriesVmMeta5stims','KmeansAccMeanVmMeta5stims',...
%     'KmeansAccSeriesVmMeta5stimsBaseline','KmeansAccMeanVmMeta5stimsBaseline'};

% this is a condensed version of KmeansUniqueness without some of the 
% extra checking and stuff i plotted and tried surrounding the analysis

getnoisefloor = 0;
binsizeVec = [1,100,400,600,800,1000:1000:10000];
% exptnames = {
%     'KP_B130_131120_p1c2c.mat'};
% numestimates = 10;
numsubsamp = 5;
numrep = 100;
numbinind = 1;
MaxTotKeyacc = [];
MaxAvgKeyacc = [];
MaxStimacc = [];
exptunqinds = [];
exptstimnames = [];
domean = 1;
getbaseline = 0;

for icell=1:size(dat,2)
    %only use one unique stimulus per cell (no doubling on dB SPL or using
    %the warped tempo stimuli
   
    loudstim = [];
    not_db_ws = [];
    notwarpstim = [];
    for istim = 1:size(stimcond,2) %grab only loudest db of db stims... or non-db stims from this set of stimuli
        %         if isempty(regexp(stimcond(istim).wavnames,'ws'))
        %             notwarpstim = [notwarpstim,istim];
        %         end
        if isempty(regexp(stimcond(istim).wavnames,'d')) && isempty(regexp(stimcond(istim).wavnames,'ws'))
            not_db_ws = [notdbstim, istim];
        end
        
        if ~isempty(regexp(stimcond(istim).wavnames,'d'))
            if ~isempty(regexp(stimcond(istim).wavnames,'d80')) || ~isempty(regexp(stimcond(istim).wavnames,'db80'))
                loudstim = [loudstim,istim];
            end
        end
    end
    
    
    
    stiminds = union(not_db_ws, loudstim); %left with only "unique" stimuli....
    %     stiminds = union(stiminds, notwarpstim);
    if isempty(stiminds)
        continue
    end
    allnames = [];
    for istim = 1:max(size(stiminds))
        thisind = stiminds(istim);
        allnames{istim} = stimcond(thisind).wavnames;
    end
    exptstimnames{iexpt} = allnames;

    %     for iestimates = 1:numestimates
    for ib = 1:size(binsizeVec,2)
        %         numsubsamp = 5;
        binsize = binsizeVec(ib);
        IDXhat = [];
        Chat = [];
        X = [];
        XmeanV = [];
        Xavg = [];
        xind = 1;
        %bin the response for that stimulus into substimuli of "binsize" long
        for istim = 1:max(size(stiminds))
            thisind = stiminds(istim);
            stimwav = stimcond(thisind).wavs; %this is the stimulus wavform in case you want to see it
            %(yes... it is a repeat of a motif, but that could end up
            %being something interesting to note... because the
            %response during the second repeat is usually different
            %from the first repeat even though the stimulus is exactly
            %the same... context/history effects
            stimexpt = filtesweeps(vmexpt,0,'wavnames',stimcond(thisind).wavnames); %filtering out all trials that is not that stimulus
            sigon = expt.analysis.params.baselinewin(2);
            sigoff = sigon+round(size(stimcond(thisind).wavs,1)/44100/expt.wc.dt);
            sigdata = medfilt1(stimexpt.wc.data,200,[],2);
            if getbaseline == 0
                sigdata = sigdata(:,sigon:sigoff); %all trials in response to that stimulus
            end
            if getbaseline == 1
                sigdata = sigdata(:,expt.analysis.params.baselinewin(1):sigon);
            end
            
            
            if getnoisefloor == 1
                for itrial = 1:size(sigdata,1)
                    sigdata(itrial,:) = sigdata(itrial,randperm(size(sigdata,2)));
                end
            end
            
            ntrials = size(sigdata,1);
            if getbaseline == 0
            numbins = floor((sigoff-sigon)/binsize); %how many subresponses of "binsize" width can there be
            end
            if getbaseline == 1
            numbins = floor((sigon-expt.analysis.params.baselinewin(1))/binsize); %how many subresponses of "binsize" width can there be
            end
            bins = [1:binsize:numbins*binsize]; %make the edges of the subresponse bins
            for ibin = 1:numbins-1 %parse the response into subresponses
                if domean == 0;
                    X = [X;sigdata(:,bins(ibin):(bins(ibin+1)-1))]; %stack all of the trials of all of the subresponses on top of each other
                end
                if domean == 1;
                    X = [X; mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)),2)];
                end
                IDXhat = [IDXhat;repmat(xind,size(sigdata,1),1)]; %keep track of what bin each trial "belongs to"
                [idx,c] = kmeans(sigdata(:,bins(ibin):(bins(ibin+1)-1)),1); %get an estimate of the "real" cluster center
                Chat (ibin,:) = c; %keep track of the centroid for each set of trials for each bin
                Xavg = [Xavg;mean(sigdata(:,bins(ibin):(bins(ibin+1)-1)))]; %this is the mean Vm response for each bin
                xind = xind +1;
            end
        end
        
        %         ik = numbins-1;
        
        replicates = 10; %how many times to run kmeans
        
        %subsample X
        
        if size(unique(IDXhat),1) < numsubsamp
            continue
        end
        
        hits = [];
        stimacc= [];
        keyacc = [];
        totalkeyacc = [];
        avgkeyacc = [];
        for iestimate = 1:numrep
            allstims = unique(IDXhat);
            subsamps = randsample(allstims,numsubsamp);
            subsampinds = [];
            for isub = 1:size(subsamps,1)
                subsampinds = [subsampinds;find(IDXhat == subsamps(isub))];
            end
            Xsub = X(subsampinds,:);
            IDXhatsub = IDXhat(subsampinds,:);
            [IDX, C, sumd,D] = kmeans(Xsub,numsubsamp,'replicate',replicates); %,'start',...
            %repmat(Ctemplate,[1,1,replicates]));
            
            allclusters = unique(IDX);
            allstims = subsamps;
            %             allstims = [1:numsubsamp];
            
            allKeys = perms(allclusters);
            
            for ikey = 1:size(allKeys,1)
                thiskey = allKeys(ikey,:);
                for istim = 1:size(allKeys,2)
                    hitKey = thiskey(istim);
                    thisstim = allstims(istim);
                    clusterinds = find(IDX == hitKey);
                    hitinds = find(IDXhatsub(clusterinds)==thisstim);
                    hits(ikey,istim,iestimate) = size(hitinds,1);
                    stimacc(ikey,istim,iestimate) = hits(ikey,istim,iestimate)/size(clusterinds,1);
                end
                totalkeyacc(ikey,iestimate) = sum(hits(ikey,:,iestimate))/size(IDX,1);
                avgkeyacc(ikey,iestimate) = mean(stimacc(ikey,:,iestimate));
            end
            
            
        end
        MaxTotKeyacc{iexpt,ib} = max(totalkeyacc);% max(mean(totalkeyacc,2));
        MaxAvgKeyacc{iexpt,ib} = max(avgkeyacc); %max(mean(avgkeyacc,2));
        %         MaxStimacc(ib) = mean(max(max(stimacc)));
    end
    %     end
%     matname = 'KmeansAccSeriesVmMeta5stims';
    matname = 'KmeansAccMeanVmMeta5stims';
%     matname = 'KmeansAccSeriesVmMeta5stimsBaseline';
%     matname = 'KmeansAccMeanVmMeta5stimsBaseline';
    if getnoisefloor == 1
        matname = 'KmeansAccVmMeta5StimsNoise';
    end
%     save([r.Dir.Expt '/Analysis/KmeansAccVm/' matname '.mat'], ...
%         'MaxTotKeyacc', 'MaxAvgKeyacc','binsizeVec')
end



%% plot for each bin, the mean accuracy across cells

figure;
set(gcf,'Color',[1,1,1])
hold on

%go through each of these data mats created for stim/baseline and
%mean/timeseries
%plot each as different color
matnames = {'KmeansAccSeriesVmMeta5stims','KmeansAccMeanVmMeta5stims',...
    'KmeansAccSeriesVmMeta5stimsBaseline','KmeansAccMeanVmMeta5stimsBaseline'};
% 'KmeansAccMeanVmMeta5stims'
% 'KmeansAccSeriesVmMeta5stimsBaseline'
% 'KmeansAccMeanVmMeta5stimsBaseline'
colors = {[1,0,1];[0,1,0];[0,0,0];[0.5, 0.5, 0.5]};%[0.5, 0.5, 0.5];%'k';%[0.8, 0.8, 0.8];%[0.3 0.3 0.3];

for imat = 1:size(matnames,2)
    load([analfolder matnames{imat} '.mat']);
    color = colors{imat};
    
    
    binmag = [];
    for ib = 1:size(MaxTotKeyacc,2)
        mag = [];
        exptind = 1;
        for iexpt = 4 %1:size(MaxTotKeyacc,1)
            if isempty(MaxTotKeyacc{iexpt,ib})
                mag(exptind) = nan;
                exptind = exptind +1;
            end
            if ~isempty(MaxTotKeyacc{iexpt,ib})
                mag(exptind) = mean(MaxTotKeyacc{iexpt,ib});
                %             scatter(binsizeVec(ib)*expt.wc.dt,mean(MaxTotKeyacc{iexpt,ib}),100,color,'fill')
                exptind = exptind +1;
            end
        end
        binmag(ib) = nanmean(mag);
    end
    scatter(binsizeVec(1:size(MaxTotKeyacc,2))*dat{1,1}{1}.dt,binmag,100,color,'fill')
end
    ylabel('clustering hit rate','FontSize',18)
    xlabel('window size (seconds)','FontSize',18)
    set(gca,'YLim',[0,1])
    title([num2str(numsubsamp) 'stimuli   ' num2str(numrep) 'estimates'], 'FontSize',14)
    

%% Membrane potential variance analysis 

%%%%%%%%%%%%%%%%%%%%%%
% aligned stimulus onset
%%%%%%%%%%%%%%%%%%%%%%

%each cell:signal pair is an average; the figures are the median across
%those variance vectors
vm_var = [];
vm_mean = [];
siglen = min(alloff - allon);
sigind = 1;
prestim_time = 1; %seconds
prestim_len = round(prestim_time/expt.wc.dt);
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    keepsigs=reprequire(table,trials);
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed(keepsigs));
    
    
    for isig = 1:size(thiscond,2)
        %skip warped shortened tempo stims
        if ~isempty(regexp(thiscond(isig).wavnames,'ws'))
            continue
        end
        sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        
        [sigon,sigoff]=GetSigTimes(sigexpt,thiscond,isig);
        startind = sigon - prestim_len;
        stopind = sigon + siglen;
        
         baselinevar = median(var(sigdata_filt(:,expt.analysis.params.baselinewin(1):sigon)));
        
        vm_var(sigind,:) = var(sigdata_filt(:,startind:stopind))./baselinevar; % percent baseline var
        vm_mean(sigind,:) = mean(sigdata_filt(:,startind:stopind)); %these estimates are not referenced to prestim var
        sigind = sigind + 1;
    end
    %%%need to do this where I do an average for each cell...
    % instead of this where every cell:signal pair is its own average
end
 %these estimates are not referenced to prestim var

%get ci on var
CIvar = [];
for isamp = 1:size(vm_var,2)
    CIvar(isamp,:) = getCDFconf (vm_var(:,isamp),85);
end

xtime = ([1:size(vm_var,2)]*expt.wc.dt)-prestim_time;
figure;
line(xtime,median(vm_var),'color','k','LineWidth',3)
line(xtime, CIvar(:,1),'color',[0.5 0.5 0.5])
line(xtime, CIvar(:,2),'color',[0.5 0.5 0.5])
axis tight
set(gca,'YLim',[-0.2,1.5]);
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean variance ; 95% CI','FontSize',14)
xlabel('seconds','FontSize',14)

zeroind = find(xtime ==0)
medianMedianVar = median(median(vm_var(:,zeroind:end)'))
CI = getCDFconf(median(vm_var(:,zeroind:end)'),85)
medianCIacrossTime = median(CIvar(zeroind:end,:))
numstimuli = size(vm_var,1)
[h,p]=kstest2(median(vm_var(:,zeroind:end)'),median(vm_var(:,1:zeroind-1)'))
baselineMedianVar = median(median(vm_var(:,1:zeroind-1)'))
figure;
line(xtime,mean(vm_mean),'color','k','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean membrane potential (mV)','FontSize',14)
xlabel('seconds','FontSize',14)
%%%what is the mean Vrest across this... (look at the offset "inhibition",
%%%does it go below Vrest? ...as far as I know it IS Vrest

%calculate when var drops below mean var during pre-stim...
%select that as start point to calculate time constant of decline of var

% prestim_vars = mean(vm_var(:,1:(prestim_len)));
% [xb, prestim_vardist] =empcdf(prestim_vars);
% prestim_ci = getCDFconf (prestim_vars,99);
% 
% for isamp = (prestim_len):size(vm_var,2)
%     thisvar = mean(vm_var(:,isamp));
%     if thisvar < (prestim_ci(1))
%         lowvar(isamp) = 1;
%     else lowvar(isamp) = 0;
%     end
% end
% 
% winlen = 5;
% for isamp = 1:size(lowvar,2) - winlen
%     meanlen(isamp) = sum(lowvar(1,isamp:isamp+winlen));
% end

% crossinds = crossing(meanlen,[],4);
% signif_drop = crossinds(min(find(crossinds>sigon)));
figure;hold on
line([1:size(vm_var,2)],mean(vm_var),'color','k')
line([1:size(vm_var,2)],median(vm_var),'color','r')
title('median var in red, mean var in black')
xonset = input('xtime onset off of median')
% % fitend = round(0.5/expt.wc.dt);
% vmdrop = mean(vm_var(:,(min(signif_drop)):end));
% % vmdrop = vmdrop(1,1:fitend) - median(vmdrop);
stopt = xonset + 0.2/expt.wc.dt;
vmdrop = median(vm_var);
vmdrop = vmdrop - median(vmdrop);
vmdrop =(vmdrop(1,xonset:stopt));


vmdrop_time = [1:size(vmdrop,2)]*expt.wc.dt;
cf = fit(vmdrop_time',vmdrop','exp1');
fitline = cf.a*exp(cf.b*vmdrop_time);

figure;hold on
line(vmdrop_time,vmdrop,'color','k')
line(vmdrop_time,fitline,'color','r')

timeconstant = -1/cf.b

%%%%%%%%%%%%% this is the figure to save and plot
figure; hold on
vmdrop = median(vm_var);
xtime = ([1:size(vmdrop,2)]*expt.wc.dt) - prestim_len*expt.wc.dt;
line(xtime,vmdrop,'color','k')
% fitline = cf.a*exp(cf.b*([1:size(vmdrop,2)-prestim_len]*expt.wc.dt));
xtime = ([1:size(fitline,2)]*expt.wc.dt) - prestim_len*expt.wc.dt + xonset*expt.wc.dt;
line(xtime,fitline+median(vmdrop),'color','r')
set(gca,'YLim',[0,1.5],'XLim',[-0.5,0.5])
title(num2str(timeconstant))
%%%%%%%%%%
%%%%%%%% baseline level looks below "1"...
%%%%%% because plotting median? ...so try using median during analysis...
%%%%%%%%%%%%%%%%% would change the % that it dropped, but better than
%%%%%%%%%%%%%%%%% baseline being below 1 since baseline is the ref

prestimvar = mean(median(vm_var(:,1:xonset)))
stimvar = mean(median(vm_var(:,(xonset+round((timeconstant*6)/expt.wc.dt)):end)))
vardec = stimvar/prestimvar

%%%%%%%%%%%%%%%%%%%%%
% aligned stimulus offset
%%%%%%%%%%%%%%%%%%%%%%

%make the var plotted as a % of baseline var
off_vm_var = [];
off_vm_mean = [];
sigind = 1;
preoffset_len = min(alloff - allon);
postoffset_len = min(alldatalen - alloff);
for iexpt=1:size(repexpts,1)
    %         if ~isempty(regexp(repexpts{iexpt},unq_expts{iunq}))
    thisexpt=repexpts{iexpt};
    load([r.Dir.Expt thisexpt])
    vmexpt=filtesweeps(expt,0,'Vm',0); %filter expt for 0 mV assuming
    table=getClampTab(expt,{'clamp',0});
    thiscond=getsubstimcond(expt.stimcond,table.sigsplayed);
    startind = alloff(sigind)-preoffset_len;
    stopind = alloff(sigind) + postoffset_len;
    
    for isig = 1:size(thiscond,2)
         if ~isempty(regexp(thiscond(isig).wavnames,'ws'))
            continue
        end
        sigexpt = filtesweeps(vmexpt,0,'wavnames',thiscond(isig).wavnames);
        sigdata = sigexpt.wc.data*1000;
        sigdata_filt = medfilt1(sigdata,200,[],2);
        
        baselinevar = median(var(sigdata_filt(:,expt.analysis.params.baselinewin(1):allon(sigind))));
        
        off_vm_var(sigind,:) = var(sigdata_filt(:,startind:stopind))./baselinevar;
        off_vm_mean(sigind,:) = mean(sigdata_filt(:,startind:stopind));
        sigind = sigind + 1;
    end
    %%%need to do this where I do an average for each cell...
    % instead of this where every cell:signal pair is its own average
end

%get ci on var
for isamp = 1:size(off_vm_var,2)
    confint(isamp,:) = getCDFconf (off_vm_var(:,isamp),85);
end

xtime = ([1:size(off_vm_var,2)]*expt.wc.dt)-preoffset_len*expt.wc.dt;
figure;
line(xtime,median(off_vm_var),'color','k','LineWidth',3)
line(xtime, confint(:,1),'color',[0.5 0.5 0.5])
line(xtime, confint(:,2),'color',[0.5 0.5 0.5])
axis tight
set(gca,'YLim',[-10,90]);
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean variance ; 95% CI','FontSize',14)
xlabel('seconds','FontSize',14)

figure;
line(xtime,median(off_vm_mean),'color','k','LineWidth',3)
axis tight
ylims = get(gca,'YLim');
line([0,0],[ylims(1),ylims(2)],'LineStyle','--','color','k');
ylabel('mean membrane potential (mV)','FontSize',14)
xlabel('seconds','FontSize',14)

%calculate when var drops below mean var during pre-stim...
%select that as start point to calculate time constant of decline of var

prestim_vars = mean(off_vm_var(:,1:(prestim_len)));
[xb, prestim_vardist] =empcdf(prestim_vars);
prestim_ci = getCDFconf (prestim_vars,99);

for isamp = (prestim_len):size(off_vm_var,2)
    thisvar = mean(off_vm_var(:,isamp));
    if thisvar < (prestim_ci(1))
        lowvar(isamp) = 1;
    else lowvar(isamp) = 0;
    end
end

winlen = 10;
for isamp = 1:size(lowvar,2) - winlen
    meanlen(isamp) = sum(lowvar(1,isamp:isamp+winlen));
end
signif_drop = crossing(meanlen,[],9);

% % fitend = round(0.5/expt.wc.dt);
% vmdrop = mean(off_vm_var(:,(min(signif_drop)):end));
% % vmdrop = vmdrop(1,1:fitend) - median(vmdrop);
% vmdrop = vmdrop - median(vmdrop);
figure;hold on
% line([1:size(off_vm_var,2)],mean(off_vm_var),'color','k')
line([1:size(off_vm_var,2)],median(off_vm_var),'color','r')
title('median var in red, mean var in black')
xonset = input('xtime onset off of median')
xoffset = input('xtime back to baseline')
yattach_val = input('bottom of begin to rise (yval)')
% % fitend = round(0.5/expt.wc.dt);
% vmdrop = mean(vm_var(:,(min(signif_drop)):end));
% % vmdrop = vmdrop(1,1:fitend) - median(vmdrop);

vmrise = median(off_vm_var);
vmrise =(vmrise(1,xonset:xoffset));
% attach = nan(1,2000);
% attach(:) = yattach_val;
% vmrise = [attach vmrise];
vmrise = vmrise - yattach_val;


vmrise_time = [1:size(vmrise,2)]*expt.wc.dt;
cf = fit(vmrise_time',vmrise','exp1');
fitline = cf.a*exp(cf.b*vmrise_time);

figure;hold on
line(vmrise_time,vmrise,'color','k')
line(vmrise_time,fitline,'color','r')


timeconstant = -1/cf.b

figure; hold on
vmrise = median(off_vm_var);
xtime = ([1:size(vmrise,2)]*expt.wc.dt) - preoffset_len*expt.wc.dt;
line(xtime,vmrise,'color','k')
xtime = ([1:size(fitline,2)]*expt.wc.dt) - preoffset_len*expt.wc.dt + xonset*expt.wc.dt;
line(xtime,fitline+yattach_val,'color','r')
axis tight
set(gca,'YLim',[0,1.5],'XLim',[-0.5,1])

%% spiking variance analysis
% use 'SpikingReproducibility' script for spiking analysis


%% effect of inhibition on spiking responses
% use 'SpikingResponse_InhibitionEffect' 



















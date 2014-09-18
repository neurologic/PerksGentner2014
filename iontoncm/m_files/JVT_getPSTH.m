% 10feb2014 - KP
%     -changed stimpath to work on local computer Krista's Macbook Pro
%     -changed binsize from 50msec to 1msec to work with FanoFactorAnalysis
% 
%   -recreated this to make it output psth as structure for FanoFactor
%   analysis to work on insted of plotting...

function [base_spikes, gz_spikes] = JVT_getPSTH(toelist, bps, start,stop, blocksize,blocktypes,sitetype, smooth, doprintout,dosave)
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};

% binsize = 1; %in msec  (this is like "dt")
% bps = 1000/binsize;
% bps = 10000;

%stimpath = '~/stims/gng_motifs';
stimpath = '/Users/kperks/GitHub/iontoncm/forstrf/strfstims';


clear toes;
clear n;
clear psth;
clear nspikes;

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
maxall=[];
%read in the spike2 format toelist
[stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);
count = mean(nspikes);

stim = strtok(stimfile, '_');
stim = strtok(stim, '.');
stimdur = stimtimes(strmatch(stim, stimnames, 'exact'));

% don't care about a raster right now
% figure;
% for itrial = 1:size(toes,2)
%     thesespikes = cell2mat(toes{itrial});
%      for ispike = 1:size(thesespikes,1)  
%         line([thesespikes(ispike) thesespikes(ispike)], [double(itrial)-.4 double(itrial)+.4],'Linewidth',.5,'Color','k');                                      
%     end     
% end

%generate histogram

stop = stimdur;
psth=[];
xbins = linspace(start, stop, (stop-start)*bps);
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

%     if isempty(newtoes)==1
%         
%         tpsth=zeros(1,length(xbins));
%     else
%         
%         tpsth=histc(newtoes, xbins);
%     end
%     psth(moo,:) = tpsth';

end

% spikes = logical(psth);
spikes = (psth);

base_spikes = spikes(basereps,:);
gz_spikes = spikes(gzreps,:);
% recpsth=mean(psth(recreps,:)./(binsize/1000));
% 
% alpha(0.8)
% 
% 
% maxall=max([basepsth gzpsth recpsth]);
% 
% if isempty(find(blocktypes==3))==0
%     gbpsth=mean(psth(gbreps,:)./(binsize/1000));
%     plot(xbins,gbpsth,'m','linewidth',1)
% end
% 
% if isempty(find(blocktypes==4))==0
%     sacpsth=mean(psth(sacreps,:)./(binsize/1000));
%     plot(xbins,sacpsth,'b','linewidth',1)
% end
% 
% if isempty(find(blocktypes==24))==0
%     gzsacpsth=mean(psth(gzsacreps,:)./(binsize/1000));
%     plot(xbins,gzsacpsth,'y','linewidth',1)    
%     maxall=max([maxall gzsacpsth]);
% end
% 
% if isempty(find(blocktypes==23))==0
%     gzgbpsth=mean(psth(gzgbreps,:)./(binsize/1000));
%     plot(xbins,gzgbpsth,'color',[.5 .5 .5],'linewidth',1);
%     
% end
% 
% 


end  %function end

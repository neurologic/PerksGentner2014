% 10feb2014 - KP
%     -changed stimpath to work on local computer Krista's Macbook Pro
%     -changed binsize from 50msec to 1msec to work with FanoFactorAnalysis


function [ax]=rasterionto(toelist, start,stop, blocksize,blocktypes,sitetype, smooth, doprintout,dosave)
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};
figure
binsize = 50; %in msec
bps = 1000/binsize;
%stimpath = '~/stims/gng_motifs';
stimpath = '/Users/kperks/GitHub/iontoncm/forstrf/strfstims';
ax(3)=0;

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

figname = strrep(sprintf('%s  %s  %s  %s  %s  %s', subjectID,pen,site,sort,stimfile,sitetype), '_', '\_');
stim = strtok(stimfile, '_');
    stim = strtok(stim, '.'); 
       stimdur = stimtimes(strmatch(stim, stimnames, 'exact'));
%stop=stimdur+2;
% =12;
% generate the raster
ax(1) = subplot(3,1,1);
hold on; 
xmin=start;xmax=stop;
ymin=0.40;
ymax=double(nreps)+.70;
%patch([xmin xmin xmax xmax],[basereps nreps-washreps nreps-washreps basereps],[.5 .5 .5],'edgealpha',0)
%patch([xmin xmin xmax xmax],[basereps basereps+drugreps basereps+drugreps basereps],[.5 .5 .5],'edgealpha',0)
for i=1:nreps
    currRep=cell2mat(toes{i});
    for j = 1:length(currRep)
        line([currRep(j) currRep(j)], [double(i)-.4 double(i)+.4],'Linewidth',.5,'Color','k');
       %    plot(cell2mat(toes{:,i}), ones(nspikes(i),1).* double(i), 'k.','MarkerSize',1);                                                 
                                                              
    end                                                            
end
hold  off;

axis([xmin xmax ymin ymax]);
ylabel('Rep');
ax_pos = get(ax(1), 'Position');
set(ax(1),'XTickLabel', '');
title(figname);
set(gca,'xtick',[]);
set(gca,'yaxislocation','left');
set(gca,'tickdir','out');
box on;

%generate histogram
ax(2) = subplot(3,1,2);
ax2_pos = get(ax(1), 'Position');
ax2_pos(2)= ax2_pos(2)-ax2_pos(4);
set(ax(2), 'Position',ax2_pos); 
nbins = (stop-start) * bps;

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
   
end

basepsth=mean(psth(basereps,:)./(binsize/1000));
gzpsth=mean(psth(gzreps,:)./(binsize/1000));
recpsth=mean(psth(recreps,:)./(binsize/1000));
hold on;
h=area(xbins,basepsth);
set(h,'facecolor','k')
alpha(0.8)
%    h=area(xbins,gzpsth);
%    set(h,'facecolor','r','edgecolor','w')
%    alpha(0.8)
%    h=area(xbins,recpsth);
%    set(h,'facecolor','g','edgecolor','w')
%    alpha(0.8)
%
plot(xbins,gzpsth,'r','linewidth',1)
plot(xbins,recpsth,'g','linewidth',1)

maxall=max([basepsth gzpsth recpsth]);

if isempty(find(blocktypes==3))==0
    gbpsth=mean(psth(gbreps,:)./(binsize/1000));
    plot(xbins,gbpsth,'m','linewidth',1)
end

if isempty(find(blocktypes==4))==0
    sacpsth=mean(psth(sacreps,:)./(binsize/1000));
    plot(xbins,sacpsth,'b','linewidth',1)
end

if isempty(find(blocktypes==24))==0
    gzsacpsth=mean(psth(gzsacreps,:)./(binsize/1000));
    plot(xbins,gzsacpsth,'y','linewidth',1)    
    maxall=max([maxall gzsacpsth]);
end

if isempty(find(blocktypes==23))==0
    gzgbpsth=mean(psth(gzgbreps,:)./(binsize/1000));
    plot(xbins,gzgbpsth,'color',[.5 .5 .5],'linewidth',1);
    
end
  
   ymax=maxall+0.1*maxall
   
  % ymax=130;
   
   axis([xmin xmax 0 ymax]);
ylabel('spikes/sec');
set(ax(2),'XTickLabel', '');
set(gca,'xtick',[]);
set(gca,'yaxislocation','left');
set(gca,'tickdir','out');
box on;
hold off   
   %figure
   

%do the spectrogram
 ax(3) = subplot(3,1,3);
 ax3_pos = get(ax(3), 'Position');
 ax3_pos(4)= ax2_pos(2)-ax3_pos(2);
 set(ax(3), 'Position',ax3_pos);

 fullstim = [stimpath '/' stimfile];
 [Y,FS,nbits]=wavread(fullstim);
 nsamples = length(Y);
%now pad the sound so that it is aligned with the raster and psth
 if (start<0)
     Y = [zeros(abs(start)* FS,1); Y];
    % fprintf('adding %g samples to the start of the stim\n', abs(start)*FS);
 elseif (xmin>0)
     Y = Y((start*FS)+1:nsamples); 
    % fprintf('subtracting %g samples from the front of the stimuli\n', start*FS);
 end

 if (stop*FS>nsamples)
     Y = [Y; zeros((stop*FS)-nsamples,1) ];
    % fprintf('adding %g samples to the end of the stimulus\n', stop*FS-nsamples);
 else
     Y = Y(1:((stop-start)*FS));   
   %  fprintf('subtracting %g samples from the end of the stimulus\n', nsamples-(stop*FS));
 end

% then do the spectrogram
 nfft = 512;
 window = 512;
 overlap = .75;
 nlap = round(window*overlap);
 [S,F,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');
FS;
 newT = xmin:((xmax-xmin)/length(T)):xmax;
 freqs=0:22050/(nfft/2):22050;
  if(length(newT)~=length(newT))
     newT=newT(1:length(T));
  end
  clim = [-200  -65];
 imagesc(newT,freqs, 20*log10(P),clim);  
 axis xy;
 colormap(flipud(colormap(gray(256))));
 axis([xmin xmax 0 11025]);
 xlabel('Time(sec)');
 ylabel('Freq (Hz)');
set(gca,'yaxislocation','left');
set(gca,'xaxislocation','bottom');
set(gca,'tickdir','out');
set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [360 635 883 286]);
set(gcf,'Name',figname,'NumberTitle','off' );
%set(gcf, 'PaperPositionMode', 'auto', 'Toolbar', 'none', 'MenuBar', 'none');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperPosition', [.15 .15 10.5 8]);

if (doprintout)
  %set(gcf,'PaperOrientation', 'landscape');
  print -dpsc2 -cmyk % to print in color
  %print -noui 
end

if dosave==1
    outfname = sprintf('%s%s%s',subjectID,site,stim);      
    saveas(gcf, outfname,'epsc');
end



end  %function end

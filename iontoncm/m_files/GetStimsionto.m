% 10feb2014 - KP
%     -changed stimpath to work on local computer Krista's Macbook Pro
%     -changed binsize from 50msec to 1msec to work with FanoFactorAnalysis


function [fullstim, Y, FS]=GetStimsionto(toelist, start,stop, blocksize,blocktypes,sitetype, smooth, doprintout,dosave)
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};

binsize = 50; %in msec
bps = 1000/binsize;
%stimpath = '~/stims/gng_motifs';
stimpath = '/Users/kperks/GitHub/iontoncm/forstrf/strfstims';
ax(3)=0;

clear toes;
clear n;
clear psth;
clear nspikes;


%read in the spike2 format toelist
[stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);

stim = strtok(stimfile, '_');
    stim = strtok(stim, '.'); 
       stimdur = stimtimes(strmatch(stim, stimnames, 'exact'));


 fullstim = [stimpath '/' stimfile];
 [Y,FS,nbits]=wavread(fullstim);

% % then do the spectrogram
%  nfft = 512;
%  window = 512;
%  overlap = .75;
%  nlap = round(window*overlap);
%  [S,F,T,P] = spectrogram(Y,window,nlap,nfft,FS, 'yaxis');
% FS;
%  newT = xmin:((xmax-xmin)/length(T)):xmax;
%  freqs=0:22050/(nfft/2):22050;
%   if(length(newT)~=length(newT))
%      newT=newT(1:length(T));
%   end
%   clim = [-200  -65];
%  imagesc(newT,freqs, 20*log10(P),clim);  
%  axis xy;
%  colormap(flipud(colormap(gray(256))));
%  axis([xmin xmax 0 11025]);
%  xlabel('Time(sec)');
%  ylabel('Freq (Hz)');
% set(gca,'yaxislocation','left');
% set(gca,'xaxislocation','bottom');
% set(gca,'tickdir','out');


end  %function end

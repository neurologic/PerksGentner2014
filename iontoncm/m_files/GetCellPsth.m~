function [Psth]=GetCellPsth(toelist, blocksize,blocktypes,sitetype)
stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940 10.571 9.932 1.000 1.000 1.000 1.000 1.000 1.000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000 5.0000];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'; 'H1'; 'H2'; '2'; '4'; '6'; '8'; '10'; '12'; 'A2short'; 'A4short'; 'B1short'; 'B2short'; 'D1short'; 'D2short'; 'E1short'; 'E2short'};
binsize = 50; %in msec
bps = 1000/binsize;

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

files = dir('*dbmean_toe.txt');
for i=1:length(files)
    toelist = files(i).name;
    %read in the spike2 format toelist
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);

    stim = strtok(stimfile, '_');
    stim = strtok(stim, '.'); 
    stimdur = stimtimes(strmatch(stim, stimnames, 'exact'));
    stop=stimdur;
    start=0;

psth=[];
xbins = linspace(start, stop, (stop-start)*bps);    
    for moo=1:nreps
            
            newtoes = toes{moo}{1}(find(toes{moo}{:}>=0));%only use spikes after stim starts
            newtoes = newtoes(find(newtoes<=stimdur)); %get rid of spikes after stim over
                        
            if isempty(newtoes)==1
                
                tpsth=zeros(1,length(xbins));
            else 
                
                tpsth=histc(newtoes, xbins);
            
                
                
            end
            psth(moo,:) = tpsth';
            
            
         
    end
  
   Psth.Base=mean(psth(basereps,:)./(binsize/1000));
   Psth.Gz=mean(psth(gzreps,:)./(binsize/1000));
   Psth.Rec=mean(psth(recreps,:)./(binsize/1000));
   
   
   
end

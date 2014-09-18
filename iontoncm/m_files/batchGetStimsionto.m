
%will make rasters and a plot of each for all stims in a site
function [wavname, wav, FS] = batchGetStimsionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
%files = dir('*_warp*toe.txt'); % to work on warped responses
%files = dir('*mot*_toe.txt'); %to work on responses to single motifs
% %files = dir('*_1_*_toe.txt'); %to work on baseline responses
% files = dir('*_toe.txt');
% for i=1:length(files)
%  if files(i).isdir==0
%   plot_rasterionto(files(i).name, -2, 10,10, 0,0);
%  end 
% end
doprintout = 0;
dosave = 0;

files=dir('*dbmean_toe.txt');

wavname = [];
wav = [];

for i=1:length(files)
    
    [stimfile,Y,FS] = GetStimsionto(files(i).name,start,stop, blocksize,blocktype,sitetype,smooth,doprintout,dosave);
    wavname{i} = stimfile;
    wav{i} = Y;
    
end



end





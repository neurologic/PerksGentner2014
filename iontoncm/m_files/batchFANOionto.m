function [base_data, gz_data] = batchFANOionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
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

for i=1:length(files)
    close all
    clearvars -except PMDdata data base_data gz_data dosave files i start stop blocksize blocktype sitetype smooth doprintout
    
%  batchrasterionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
%  [base, gz]= FanoFormat_ionto(files(i).name,-2,12,10,[1 2 5],'mu',0,0);
 [base, gz]= FanoFormat_ionto(files(i).name,start,stop,blocksize,blocktype,sitetype,smooth,doprintout);
 
 base_data(i).spikes = base;
 gz_data(i).spikes = gz;

end
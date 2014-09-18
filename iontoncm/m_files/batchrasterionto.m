
%will make rasters and a plot of each for all stims in a site
function batchrasterionto(start,stop,blocksize,blocktype,sitetype,smooth,doprintout)
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
    clearvars -except dosave files i start stop blocksize blocktype sitetype smooth doprintout
    
    
    
    rasterionto(files(i).name,start,stop, blocksize,blocktype,sitetype,smooth,doprintout,dosave)
%     psth = FanoFormat_ionto(files(i).name,start,stop, blocksize,blocktype,sitetype,smooth,doprintout,dosave)
    
    
end


% clear all; close all; plot_rasterionto2('ss001A2_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001A4_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001B1_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001B2_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001D1_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001D2_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001E1_68dbmean_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss001E2_68dbmean_toe.txt',-2,5,0,1)
% 
% clear all; close all; plot_rasterionto2('ss00110_khz_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss00112_khz_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss0012_khz_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss0014_khz_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss0016_khz_toe.txt',-2,5,0,1)
% clear all; close all; plot_rasterionto2('ss0018_khz_toe.txt',-2,5,0,1)
% 

end





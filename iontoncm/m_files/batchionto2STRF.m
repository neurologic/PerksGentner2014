!rm _STRT.txt

files=dir('*stim.mat');
for n=1:length(files)
    ionto2STRFfile(files(n).name)
    
end
    
    

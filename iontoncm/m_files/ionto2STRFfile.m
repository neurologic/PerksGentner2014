%TOE2STRFFILE.M
%
%Convert ToeFiles to STRFPak-readable response file format
%Overwrites any existing <filename>_STRF.txt file already created with this
%function.
%
%Input: 
%  toefile: path to toefile to convert
%
%Output:
%  responseFileName: filename of converted response file. This file is
%  suitable for using with STRFPak.


%JMJ 08/06/07
function responseFileName = ionto2STRFfile(file)
srfactor=44100/43868;
load(file)

tip=strtok(file,'.');

responseFileName = [tip '_STRF.txt'];

%delete(responseFileName);

for i = 1:length(spikes),
    dlmwrite(responseFileName,spikes{i}'*srfactor*1000,'delimiter', ' ','-append');
end
end
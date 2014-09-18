function [motifsByRep,motifsByRep_shuffled, minReps]=GetMotFRsJTionto(drugclass)


%codes for motif identity for each mot in each song
mA1 = [1 1 1 2 2 2 2 3 3 3 3 ];
mA2 = [4 4 5 5 5 5 5 6 6 6 6 7];
mA3 = [8 9 9 9 10 10 10 10 10 11 11 3 3];
mA4 = [48 49 49 4 4 4 5 5 5 5 5 52];
mA5 = [53 53 53 54 54 54 55 55 55 55];
mB1 = [13 13 13 14 14 14 15 16 16 17 17 18 18];
mB2 = [19 19 19 20 20 21 21 21 21 22 22 22 23];
mB3 = [24 25 26 26 27 28];
mB4 = [56 57 57 58 59 59 59 60 60 60 61 61];
mB5 = [62 63 64 13 13 13 14 14];
mC1 = [32 33 33 33 34 34 35 35 36 37 38 38 39];
mC2 = [40 41 41 41 42 43 43 44 45 45 46 47];
mD1 = [65 66 67 68 68 69 70 70 71 72 73 74 75];
mD2 = [76 77 78 79 79 80 80 80 80 80 81 82 83 83];
mE1 = [84 85 86 87 87 88 89 90 91 92 93 93 93 94 95];
mE2 = [96 96 97 98 99 92 93 93 100 94 101 101 102];
mN4 = [103 104 105 105 105 106 107 108 108 109 109 110 110 111];
mN6 = [118 119 120 121 121 122 123 123 124 124 125];
mN3 = [126 127 128 129 130 131 132 133 134 135 136 137 138];
mN2 = [32 33 33 33 34 34 35 35 36 37 38 38 39];
mN1 = [139 140 141 142 143 144 145 145 146 147 147 148 149];
mN5 = [40 41 41 41 42 43 43 44 45 45 46 47];

%break points b/t each motif in each song
cA1 = [0 0.809 1.668 2.567 3.582 4.544 5.533 6.603 7.461 8.343 9.166 9.914];
cA4 = [0 .622 1.531 2.401 3.300  4.198  4.952  5.864  6.704 7.394  8.102  8.801 9.215];
cA2 = [0  .78  1.59  2.45  3.25  4.06  4.81  5.58  6.32 7.01  7.69  8.50 9.0756];
cA3 = [0 0.711 1.724 2.525 3.420 4.261 5.123 5.995 6.680 7.360 8.020 8.818 9.675 10.432];
cA5 = [0 1.154 2.474 3.432 4.214 5.074 5.942 6.632 7.504 8.330 9.272];
cB1 = [0 .731  1.698  2.707  3.456  4.228  5.247  5.996  6.536  7.030  7.614  8.300 8.916 9.571]; 
cB2 = [0  0.906  1.594  2.449  3.220  4.107  4.827  5.539  6.296  6.884  7.437  8.010 8.678];
cB3 = [0 1.129 3.079 4.865 6.679 8.153 9.290];
cB4 = [0 0.665 1.615 2.867 3.944 4.588 5.352 6.173 6.951 7.714 8.385 9.176 9.652];
cB5 = [0 1.481 2.730 4.706 5.578 6.584 7.640 8.418 9.076];
cC1 = [0  0.751  1.570  2.335  3.152  3.800  4.540  5.295  6.104 6.999  8.011  8.577  9.361 10.523];
cC2 = [0  0.731  1.325  2.140  3.012  3.510  4.787  6.228  7.153  7.939  8.785  9.779 10.707];
cD1 = [0 0.844 1.749 2.921 4.312 5.578 5.971 6.710 7.423 8.185 8.481 9.140 9.491 10.122];
cD2 = [0 0.787 1.333 2.081 2.907 3.838 4.712 5.65 6.305 6.80 7.558 8.142 8.681 9.475 10.151];
cE1 = [0 0.637 1.440 2.193 2.739 3.605 4.042 4.543 5.087 5.871 6.547 7.443 8.208 8.903 9.350 9.950];
cE2 = [0 1.138 2.97 3.440 4.177 4.904 5.694 6.654 7.347 7.943 8.570 9.320 10.150];
cN1 = [0 1.094 2.129 3.026 3.745 4.377 5.117 5.889 6.799 7.339 8.044 8.678 9.317 10.251];
cN2 = [0 0.828 1.646 2.380 3.194 3.676 4.571 5.366 6.156 7.041 8.013 8.576 9.364 10.523];
cN3 = [0 1.550 2.563 3.271 4.337 5.223 6.006 6.454 7.192 7.954 8.483 8.993 9.491 10.136];
cN4 = [0 0.405 1.375 1.974 2.555 3.130 3.563 4.719 5.272 5.785 6.521 7.477 8.125 9.284 10.606];
cN5 = [0 0.742 1.341 2.172 3.053 4.784 6.227 7.151 7.943 8.784 9.778 10.717];
cN6 = [0 0.971 1.853 2.977 4.321 5.576 5.981 6.722 7.434 8.200 9.119 10.122];
cS1 = [0 1.5 3.5 5.5 7.5];
cS2 = [0 3.0 7.0 10.0];

stimtimes = [ 9.915 9.082 10.432 9.225 9.273 9.579 9.334 9.292 9.654 9.078 10.257 10.532 10.135 10.611 10.719 10.124 10.532 10.719 10.124 10.160 9.954 10.155 9.000 10.940];
stimnames = {'A1'; 'A2'; 'A3'; 'A4'; 'A5'; 'B1'; 'B2'; 'B3'; 'B4'; 'B5'; 'N1'; 'N2'; 'N3'; 'N4'; 'N5'; 'N6'; 'C1'; 'C2'; 'D1'; 'D2'; 'E1'; 'E2'; 'S1'; 'S2'};

bird = {'st444'; 'st445'; 'st493'; 'st484'; 'st504'; 'st450'; 'st490'; 'st466'};
trainset = {'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'; 'E1' 'E2' 'A2' 'A4'; 'E1' 'E2' 'A2' 'A4'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'B1' 'B2'; 'D1' 'D2' 'E1' 'E2'};
goset = {'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'; 'E1' 'E2'; 'E1' 'E2'; 'B1' 'B2'; 'A2' 'A4'; 'D1' 'D2'};
nogoset = { 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'; 'A2' 'A4'; 'A2' 'A4'; 'D1' 'D2'; 'B1' 'B2'; 'E1' 'E2'};
novset = { 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'; 'B1' 'B2' 'D1' 'D2'; 'B1' 'B2' 'D1' 'D2'; 'A2' 'A4' 'E1' 'E2'; 'D1' 'D2' 'E1' 'E2'; 'A2' 'A4' 'B1' 'B2'};

  
  if strmatch(drugclass,'base')==1
      files = dir('*00*base.txt');
  elseif strmatch(drugclass,'gz')==1
    files = dir('*00*gz.txt');  
  elseif strmatch(drugclass,'rec')==1
    files = dir('*00*rec.txt');  
  end
  
minreps=100;
 crap=[];
for u=1:length(files)
    toelist = files(u).name;
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);
    if nreps<minreps
        minreps=nreps;
    end
    [t,k]=strtok(files(u).name,'_');
    [y,z]=strtok(t,sort);
    if strmatch(y,'S1')==1
       crap=[crap u];
    elseif strmatch(y,'S2')==1
        crap=[crap u];
    end
end
 
files(crap)=[];
subindx=strcmpi(bird,subjectID);
 %

trainstims=[];
novstims=[];
gostims=[];
nogostims=[];

allnmotrate=[];
alltmotrate=[];
allgomotrate=[];
allnogomotrate=[];
 x=1;
for u=1:length(files)
  
    toelist = files(u).name;
    [stimfile, subjectID, pen, site, sort, nreps, nspikes, toes, alltoes] = readtoe(toelist);
    stim = strtok(stimfile, '_');
       stim = strtok(stim, '.'); 
    if sum(strcmpi(stim,trainset(subindx,:)))>0
        trainstims=[trainstims  u];
    end
    if sum(strcmpi(stim,novset(subindx,:)))>0
        novstims=[novstims  u];
    end
        if sum(strcmpi(stim,goset(subindx,:)))>0
        gostims=[gostims  u];
    end
    if sum(strcmpi(stim,nogoset(subindx,:)))>0
        nogostims=[nogostims  u];
    end
    cuts = eval(['c' stim]);
    motifid = eval(['m' stim]);
          
    nreps=minreps;
    for y=1:size(cuts,2)-1;% cycle through each motif in song
         
               
        for moo=1:nreps
        
            notoes = toes{moo}{1}(find(toes{moo}{:}>=cuts(y)));%only use spikes after stim starts
            notoes = notoes(find(notoes<=cuts(y+1))); %get rid of spikes after stim over
            morate(moo)=length(notoes) / ( cuts(y+1)- cuts(y) );
        end
          
             
             
             
              if sum(strcmpi(stim,trainset(subindx,:)))>0
                    
                    alltmotrate=[alltmotrate morate'];
              end
                if sum(strcmpi(stim,novset(subindx,:)))>0
                    
                    allnmotrate=[allnmotrate morate'];
                end
                 if sum(strcmpi(stim,goset(subindx,:)))>0
                    
                    allgomotrate=[allgomotrate morate'];
                 end
                 if sum(strcmpi(stim,nogoset(subindx,:)))>0
                    
                    allnogomotrate=[allnogomotrate morate'];
                end
                
                
                
             x=x+1;
             
        
    end
end

minReps=double(minreps);

motifsByRep.go = allgomotrate;
motifsByRep.nogo = allnogomotrate;
motifsByRep.novel = allnmotrate;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get Shuffled Motif Firing Rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    class0MotifsByRepShuffled = zeros(minReps, size(allnmotrate,2));
    class1MotifsByRepShuffled = zeros(minReps, size(allgomotrate,2));
    class2MotifsByRepShuffled = zeros(minReps, size(allnogomotrate,2));

    repsByMotif = motifsByRep.go;
    for t = 1:size(repsByMotif,1),
        ordering = randperm(size(repsByMotif,2));
        class1MotifsByRepShuffled(t,:) = repsByMotif(t,ordering);
    end
    repsByMotif = motifsByRep.nogo;
    for t = 1:size(repsByMotif,1),
        ordering = randperm(size(repsByMotif,2));
        class2MotifsByRepShuffled(t,:) = repsByMotif(t,ordering);
    end
    repsByMotif = motifsByRep.novel;
    for t = 1:size(repsByMotif,1),
        ordering = randperm(size(repsByMotif,2));
        class0MotifsByRepShuffled(t,:) = repsByMotif(t,ordering);
    end

    motifsByRep_shuffled.go = class1MotifsByRepShuffled;
    motifsByRep_shuffled.nogo = class2MotifsByRepShuffled;
    motifsByRep_shuffled.novel = class0MotifsByRepShuffled;
    
    
    
end
 












        

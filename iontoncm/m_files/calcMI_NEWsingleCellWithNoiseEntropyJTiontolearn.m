%This is used to calculate MI, entropy, and noise entropy for single cells
%for Fig 2 in Nature submission. Lives in
%Code/Information/ComputationScripts/ directory.

%Calculates MI between firing rate and motif ID within each class (go,
%nogo, novel).
%jmj 3/11/09

%%


NCMMIcellIDs= { '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z1230/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z1460/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z1590/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z1810/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z2570/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z2680/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-300_L1000/Z3200/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-500_L1000/Z2000/ss001' ; '/home/jason/data/iontoncm/B444/Left_R-500_L1000/Z2450/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z1550/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z1670/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z1750/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z1860/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z2150/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z2330/ss001' ; '/home/jason/data/iontoncm/B445/Left_R-200_L1000/Z2770/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1000/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1055/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1151/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1260/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1330/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1440/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1480/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1570/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z1700/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z2650/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z920/ss001' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z950/ss003' ; '/home/jason/data/iontoncm/B493/Right_R-350_L1100/Z950/ss004' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1110/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1150/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1260/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1310/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1440/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1520/ss003' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1640/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1760/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z2190/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z2250/ss001' ; '/home/jason/data/iontoncm/B484/Left_R-300_L1000/Z1370/ss001' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z830/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2040/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2100/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2250/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2350/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2480/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2600/ss001/' ; '/home/jason/data/iontoncm/B450/Right_R-300_L1000/Z2800/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z780/ss002/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z2150/ss002/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z2210/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z2250/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z2311/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-200_L1000/Z2430/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-350_L1000/Z2800/ss001/' ; '/home/jason/data/iontoncm/B466/Left_R-350_L1000/Z2880/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z660/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z660/ss002/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2150/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2280/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2430/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2650/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2810/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2910/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z2990/ss001/' ; '/home/jason/data/iontoncm/B504/LEFT_R-300_L1000/Z3032/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z970/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z1170/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z1461/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z1570/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z1850/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z2101/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z2220/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z2620/ss001/' ; '/home/jason/data/iontoncm/B490/Right_R-300_L1000/Z2811/ss001/'};
drugclass='gz';
numFRbins = 7;
%cellIDs = 10;

c=1;
for m = 1:length(NCMMIcellIDs),
    
    currentDirectory = NCMMIcellIDs{m};
    cd(currentDirectory);
    
    class1=1;
    class2=2;
    class0=3;
    
    iSR1A = {};
    iSR2A = {};
    iSR0A = {};
    Ent1A = {};
    Ent2A = {};
    Ent0A = {};
       
    NoiseEnt1A = {};
    NoiseEnt2A = {};
    NoiseEnt0A = {};
    
    
    
    iSR1AnovelRange = {};
    iSR2AnovelRange = {};
    iSR0AnovelRange = {};
    Ent1AnovelRange = {};
    Ent2AnovelRange = {};
    Ent0AnovelRange = {};
 
    
    
    iSR1Adummy = {};
    iSR2Adummy = {};
    iSR0Adummy = {};
    Ent1Adummy = {};
    Ent2Adummy = {};
    Ent0Adummy = {};
    
    [motifsByRepA, motifsByRep_shuffledA, minRepsA] = GetMotFRsJTiontolearn(drugclass);
    
    minFRA = min(min([motifsByRepA.go motifsByRepA.nogo motifsByRepA.novel]));
    maxFRA = max(max([motifsByRepA.go motifsByRepA.nogo motifsByRepA.novel]));
   
    binsA = linspace(minFRA, maxFRA, numFRbins);
   
    binsNovelRange = linspace(min(min(motifsByRepA.novel)), max(max(motifsByRepA.novel)), numFRbins);
    
    
    %Set number of reps to use for extrapolation:
    extrapRepsMin = minRepsA - min(floor(0.4*minRepsA), 8);
    
    
    extrapReps = extrapRepsMin:minRepsA;

    for k = 1:length(extrapReps),
        C = nchoosek(minRepsA, extrapReps(k));
        maxExtrapReps = 500;
        if C <= maxExtrapReps
            %Do all possible combinations
            combsToUse = int32(1:C);            
        else
            %Get 500 of the possible combinations
            combsToUse = int32((1:500)*C/500);
        end
        
            
        for y = 1:length(combsToUse),
            randomComb = onecomb(minRepsA, extrapReps(k), combsToUse(y));
           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Build conditional distributions of ACTUAL data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            condProbA.go = ConditionalDistributionFromReps(motifsByRepA.go(1:minRepsA,:),randomComb,binsA);
            condProbA.nogo = ConditionalDistributionFromReps(motifsByRepA.nogo(1:minRepsA,:),randomComb,binsA);
            condProbA.novel = ConditionalDistributionFromReps(motifsByRepA.novel(1:minRepsA,:),randomComb,binsA);

            %need the nan2zeros because some rows of the cond prob
            %distribution will sum to zero.
            condProbAnovelRange.go = nan2zero(ConditionalDistributionFromReps(motifsByRepA.go(1:minRepsA,:),randomComb,binsNovelRange));
            condProbAnovelRange.nogo = nan2zero(ConditionalDistributionFromReps(motifsByRepA.nogo(1:minRepsA,:),randomComb,binsNovelRange));
            condProbAnovelRange.novel = nan2zero(ConditionalDistributionFromReps(motifsByRepA.novel(1:minRepsA,:),randomComb,binsNovelRange));
       
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Build conditional distributions of DUMMY data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            condProbDummyA.go = ConditionalDistributionFromReps(motifsByRep_shuffledA.go(1:minRepsA,:),randomComb,binsA);
            condProbDummyA.nogo = ConditionalDistributionFromReps(motifsByRep_shuffledA.nogo(1:minRepsA,:),randomComb,binsA);
            condProbDummyA.novel = ConditionalDistributionFromReps(motifsByRep_shuffledA.novel(1:minRepsA,:),randomComb,binsA);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calc Joint Probability distribution of 2d data:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            jointDistA.go = condProbA.go./sum(sum(sum(condProbA.go)));
            jointDistA.nogo = condProbA.nogo./sum(sum(sum(condProbA.nogo)));
            jointDistA.novel = condProbA.novel./sum(sum(sum(condProbA.novel)));
              
            jointDistAnovelRange.go = condProbA.go./sum(sum(sum(condProbAnovelRange.go)));
            jointDistAnovelRange.nogo = condProbA.nogo./sum(sum(sum(condProbAnovelRange.nogo)));
            jointDistAnovelRange.novel = condProbA.novel./sum(sum(sum(condProbAnovelRange.novel)));
               
            
            %jointDistAdummy.go = CalcJointDistFromCondDistsTemp(condProbA.go, condProbDummyA.go, length(binsA)-1);
            %jointDistAdummy.nogo = CalcJointDistFromCondDistsTemp(condProbA.nogo, condProbDummyA.nogo, length(binsA)-1);
            %jointDistAdummy.novel = CalcJointDistFromCondDistsTemp(condProbA.novel, condProbDummyA.novel, length(binsA)-1);
            
            [iSR1A{k}(y), Ent1A{k}(y), NoiseEnt1A{k}(y)] = calcMIandEntandNoiseEnt(jointDistA.go);
            [iSR2A{k}(y), Ent2A{k}(y), NoiseEnt2A{k}(y)] = calcMIandEntandNoiseEnt(jointDistA.nogo);
            [iSR0A{k}(y), Ent0A{k}(y), NoiseEnt0A{k}(y)] = calcMIandEntandNoiseEnt(jointDistA.novel);
     
            [iSR1AnovelRange{k}(y), Ent1AnovelRange{k}(y)] = calcMIandEnt(jointDistAnovelRange.go);
            [iSR2AnovelRange{k}(y), Ent2AnovelRange{k}(y)] = calcMIandEnt(jointDistAnovelRange.nogo);
            [iSR0AnovelRange{k}(y), Ent0AnovelRange{k}(y)] = calcMIandEnt(jointDistAnovelRange.novel);
     
            
           % [iSR1Adummy{k}(y), Ent1Adummy{k}(y)] = calcMIandEnt(jointDistAdummy.go);
           % [iSR2Adummy{k}(y), Ent2Adummy{k}(y)] = calcMIandEnt(jointDistAdummy.nogo);
           % [iSR0Adummy{k}(y), Ent0Adummy{k}(y)] = calcMIandEnt(jointDistAdummy.novel);


            
            
        end
    end

    [extrapA.go(m),extrapExtrasA.go(m)] = extrapolateInformation(iSR1A, extrapReps);
    [extrapA.nogo(m),extrapExtrasA.nogo(m)] = extrapolateInformation(iSR2A, extrapReps);
    [extrapA.novel(m),extrapExtrasA.novel(m)] = extrapolateInformation(iSR0A, extrapReps);

    extrapAAll(m,class1) = extrapA.go(m);
    extrapAAll(m,class2) = extrapA.nogo(m);
    extrapAAll(m,class0) = extrapA.novel(m);
    
    
    
    [extrapAent.go(m),extrapExtrasAent.go(m)] = extrapolateInformation(Ent1A, extrapReps);
    [extrapAent.nogo(m),extrapExtrasAent.nogo(m)] = extrapolateInformation(Ent2A, extrapReps);
    [extrapAent.novel(m),extrapExtrasAent.novel(m)] = extrapolateInformation(Ent0A, extrapReps);

    [extrapAnoiseEnt.go(m),extrapExtrasAnoiseEnt.go(m)] = extrapolateInformation(NoiseEnt1A, extrapReps);
    [extrapAnoiseEnt.nogo(m),extrapExtrasAnoiseEnt.nogo(m)] = extrapolateInformation(NoiseEnt2A, extrapReps);
    [extrapAnoiseEnt.novel(m),extrapExtrasAnoiseEnt.novel(m)] = extrapolateInformation(NoiseEnt0A, extrapReps);

    
    
    [extrapAnovelRange.go(m),extrapExtrasAnovelRange.go(m)] = extrapolateInformation(iSR1AnovelRange, extrapReps);
    [extrapAnovelRange.nogo(m),extrapExtrasAnovelRange.nogo(m)] = extrapolateInformation(iSR2AnovelRange, extrapReps);
    [extrapAnovelRange.novel(m),extrapExtrasAnovelRange.novel(m)] = extrapolateInformation(iSR0AnovelRange, extrapReps);

    [extrapAentNovelRange.go(m),extrapExtrasAentNovelRange.go(m)] = extrapolateInformation(Ent1AnovelRange, extrapReps);
    [extrapAentNovelRange.nogo(m),extrapExtrasAentNovelRange.nogo(m)] = extrapolateInformation(Ent2AnovelRange, extrapReps);
    [extrapAentNovelRange.novel(m),extrapExtrasAentNovelRange.novel(m)] = extrapolateInformation(Ent0AnovelRange, extrapReps);
    
    
    
    %[extrapAdummy.go(m),extrapExtrasAdummy.go(m)] = extrapolateInformation(iSR1Adummy, extrapReps);
    %[extrapAdummy.nogo(m),extrapExtrasAdummy.nogo(m)] = extrapolateInformation(iSR2Adummy, extrapReps);
    %[extrapAdummy.novel(m),extrapExtrasAdummy.novel(m)] = extrapolateInformation(iSR0Adummy, extrapReps);

end
%save /home/jamie/lintu/Jamie/Code/Information/SavedComputations/infoNewSingleCellFinal6bin_10plus_CLM_ByStimulus extrapA extrapAAll extrapExtrasA extrapAent extrapExtrasAent extrapAnoiseEnt extrapExtrasAnoiseEnt extrapAnovelRange extrapExtrasAnovelRange extrapAentNovelRange extrapExtrasAnovelRange extrapAdummy extrapExtrasAdummy
%runANOVAstats(extrapAentNovelRange.go, extrapAentNovelRange.nogo, extrapAentNovelRange.novel)
%save infoNewSingleCellDeterministicEnt extrapAent extrapAentNovelRange extrapExtrasAent extrapExtrasAentNovelRange
%runANOVAstats(extrapAnovelRange.go, extrapAnovelRange.nogo, extrapAnovelRange.novel)




%%

figure
infomat=[extrapA.go' extrapA.nogo' extrapA.novel'];
bar(mean(infomat),'k');
hold on;
errorbar(mean(infomat),sem(infomat),'.k')
ylabel('Mean M.I. \pm S.E.M.');
set(gca,'xticklabel',{'Go'; 'NoGo';'Novel'});
title(drugclass)

figure
entmat=[extrapAent.go' extrapAent.nogo' extrapAent.novel'];
bar(mean(entmat),'k');
hold on;
errorbar(mean(entmat),sem(entmat),'.k')
ylabel('Response Entropy');
set(gca,'xticklabel',{'Go'; 'NoGo';'Novel'});
title(drugclass)

figure
noiseentmat=[extrapAnoiseEnt.go' extrapAnoiseEnt.nogo' extrapAnoiseEnt.novel'];
bar(mean(noiseentmat),'k');
hold on;
errorbar(mean(noiseentmat),sem(noiseentmat),'.k')
ylabel('Noise Entropy');
set(gca,'xticklabel',{'Go'; 'NoGo';'Novel'});
title(drugclass)

%%
values1 = [extrapAentNovelRange.go; extrapAentNovelRange.nogo; extrapAentNovelRange.novel];
values2 = [extrapAent.go; extrapAent.nogo; extrapAent.novel];

run2wayANOVAstats(values1, values2)


plotNGN(extrapA, [0 0.6], 'Mean M.I. \pm S.E.M.')
set(gca, 'YTick', [0 .1 .2 .3 .4 .5 .6]);


plotNGN(extrapAent, [0 1.7], 'Response Entropy');

figure
plotNGN(extrapAnoiseEnt, [0 1.7], 'Noise Entropy');

errorbar([1 2 3], [mean(extrapAent.go) mean(extrapAent.nogo) mean(extrapAent.novel)], [std(extrapAent.go) std(extrapAent.nogo) std(extrapAent.novel)]/sqrt(78), 'sk', 'MarkerFaceColor', 'k');
hold on
errorbar([1 2 3], [mean(extrapAnoiseEnt.go) mean(extrapAnoiseEnt.nogo) mean(extrapAnoiseEnt.novel)], [std(extrapAnoiseEnt.go) std(extrapAnoiseEnt.nogo) std(extrapAnoiseEnt.novel)]/sqrt(78), 'ok', 'MarkerFaceColor', 'k');
ylim([0.8 1.7]);
xlim([0.5 3.5]);
set(gca, 'FontSize', 11);
legend('Total Entropy', 'Noise Entropy');
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Reinforced', 'Non-Reinforced', 'Novel'});
box off
ylabel('Mean entropy \pm SEM (bits)');
set(gca, 'YTick', [0.8 1 1.2 1.4 1.6]);
set(gca, 'YTickLabel', {'0', '1.0', '1.2', '1.4', '1.6'});
rotateticklabel(gca, -45);
get(gca, 'Position');
set(gcf, 'Position', [0 0 250 400]);
set(gca, 'Position', [.2657 .21 .50 .715]);
fixTickMarks(gca);

%Plot entropy histograms
subplot(2,1,1)
difference_10 = extrapAent.go - extrapAent.novel;
range = -1:0.125:1.5;
[N,bin] = histc(difference_10,range);
h=bar(range,N,1,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'EdgeColor', [.5 .5 .5]);
hold on
plot([mean(difference_10) mean(difference_10)], [0 20], 'k')
plot([0 0], [0 20], '--k');
xlim([-1.5 1.5]);
set(gca, 'XTick', [-1.5 -0.75 0 0.75 1.5]);
set(gca, 'YTick', [0 10 20]);
set(gca, 'FontSize', 11);
set(gca, 'Layer', 'top');
%set(gca, 'XTickLabel', []);
box off
ylabel('No. of Neurons');
set(gcf, 'Position', [0 0 320 250])
fixTickMarks(gca);


subplot(2,1,2)
difference_20 = extrapAent.nogo - extrapAent.novel;
range = -1:0.125:1.5;
[N,bin] = histc(difference_20,range);
h=bar(range,N,1,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'EdgeColor', [.5 .5 .5]);
hold on
plot([mean(difference_20) mean(difference_20)], [0 20], 'k')
plot([0 0], [0 20], '--k');
xlim([-1.5 1.5]);
set(gca, 'XTick', [-1.5 -0.75 0 0.75 1.5]);
set(gca, 'YTick', [0 10 20]);
set(gca, 'FontSize', 11);
set(gca, 'Layer', 'top');
%set(gca, 'XTickLabel', []);
box off
ylabel('No. of Neurons');
xlabel('\Delta H_{Response}');
fixTickMarks(gca);


figure
h = bar([1 2], [mean(difference_10) mean(difference_20)], .6);
hold on
h2 = errorbar([1 2], [mean(difference_10) mean(difference_20)], [std(difference_10) std(difference_20)]/sqrt(96), '.k', 'MarkerSize', 1);
set(h, 'FaceColor', [.8 .8 .8]);
ylabel('Mean \DeltaH_{Response} \pm S.E.M. (bits)');
xlim([.5 2.5]);
set(gca, 'FontSize', 11);
set(gca, 'XTick', [1 2]);
set(gca, 'YTick', [0 .05 .1 .15 .2]);
set(gca, 'XTickLabel', {'Go-Novel', 'Nogo-Novel'});
box off
fixTickMarks(gca);
set(gcf, 'Position', [0 0 150 250]);



subplot(2,1,2)
[N,bin] = histc(difference_20,range);
h=bar(range,N,1.05,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'LineStyle', 'none');
hold on
plot([mean(difference_20) mean(difference_20)], [0 20], 'k')
plot([0 0], [0 20], '--k');
xlim([-1.5 1.5]);
set(gca, 'XTick', [-1.5 0 1.5]);
set(gca, 'YTick', [0 20]);
set(gca, 'FontSize', 11);
%set(gca, 'XTickLabel', []);
box off
xlabel('\Delta RS')
ylabel('No. of Neurons');



%Plot MI histograms
subplot(3,1,1)
difference_10 = extrapA.go - extrapA.nogo;
range = -.8:0.08:0.8;
[N,bin] = histc(difference_10,range);
h=bar(range,N,1,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'EdgeColor', [.5 .5 .5]);
hold on
plot([mean(difference_10) mean(difference_10)], [0 23], 'k')
plot([0 0], [0 23], '--k');
xlim([-0.8 0.8]);
ylim([0 23]);
set(gca, 'XTick', [-.8 0 0.8]);
set(gca, 'YTick', [0 10 20]);
set(gca, 'FontSize', 11);
set(gca, 'Layer', 'top');
%set(gca, 'XTickLabel', []);
box off
%ylabel('No. of Neurons');
set(gcf, 'Position', [0 0 250 250])
fixTickMarks(gca);


subplot(3,1,2)
difference_10 = extrapA.go - extrapA.novel;
range = -.8:0.08:0.8;
[N,bin] = histc(difference_10,range);
h=bar(range,N,1,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'EdgeColor', [.5 .5 .5]);
hold on
plot([mean(difference_10) mean(difference_10)], [0 23], 'k')
plot([0 0], [0 23], '--k');
xlim([-0.8 0.8]);
ylim([0 23]);
set(gca, 'XTick', [-.8 0 0.8]);
set(gca, 'YTick', [0 10 20]);
set(gca, 'FontSize', 11);
set(gca, 'Layer', 'top');
%set(gca, 'XTickLabel', []);
box off
ylabel('No. of Neurons');
fixTickMarks(gca);


subplot(3,1,3)
difference_20 = extrapA.nogo - extrapA.novel;
[N,bin] = histc(difference_20,range);
h=bar(range,N,1,'histc');
set(h, 'FaceColor', [.5 .5 .5]);
set(h, 'EdgeColor', [.5 .5 .5]);
hold on
plot([mean(difference_20) mean(difference_20)], [0 23], 'k')
plot([0 0], [0 23], '--k');
xlim([-0.8 0.8]);
ylim([0 23]);
set(gca, 'XTick', [-.8 0 .8]);
set(gca, 'YTick', [0 10 20]);
set(gca, 'FontSize', 11);
set(gca, 'Layer', 'top');
%set(gca, 'XTickLabel', []);
box off
%ylabel('No. of Neurons');
xlabel('\Delta MI (bits)');
fixTickMarks(gca);


figure
h = bar([1 2], [mean(difference_10) mean(difference_20)], .6);
hold on
h2 = errorbar([1 2], [mean(difference_10) mean(difference_20)], [std(difference_10) std(difference_20)]/sqrt(96), '.k', 'MarkerSize', 1);
set(h, 'FaceColor', [.8 .8 .8]);
ylabel('Mean \DeltaMI \pm S.E.M. (bits)');
xlim([.5 2.5]);
set(gca, 'FontSize', 11);
set(gca, 'XTick', [1 2]);
set(gca, 'YLim', [0 0.15]);
set(gca, 'YTick', [0 .05 .1 .15]);
set(gca, 'XTickLabel', {'Go-Novel', 'Nogo-Novel'});
box off
fixTickMarks(gca);
set(gcf, 'Position', [0 0 150 250]);



plotNGN(extrapA.go, extrapA.nogo, extrapA.novel, 'Information (bits)');
box off
set(gcf, 'Position', [0 0 150 400])
runANOVAstats(extrapA.go, extrapA.nogo, extrapA.novel)

plotNGN(extrapAent.go, extrapAent.nogo, extrapAent.novel, 'H_{Response} (bits)');
box off
set(gcf, 'Position', [0 0 150 400])
runANOVAstats(extrapAent.go, extrapAent.nogo, extrapAent.novel)





%Test for normality of entropy data:
[h,p] = lillietest(extrapAent.go - extrapAent.novel)
[h,p] = lillietest(extrapAent.nogo - extrapAent.novel)
[h,p] = lillietest(extrapAent.nogo - extrapAent.go)


[h,p] = lillietest(extrapA.go - extrapA.novel)
[h,p] = lillietest(extrapA.nogo - extrapA.novel)
[h,p] = lillietest(extrapA.nogo - extrapA.go)


[p] = friedman([extrapA.go' extrapA.nogo' extrapA.novel'])


iSR1(find(iSR1 < 0)) = 0;
iSR2(find(iSR2 < 0)) = 0;
iSR0(find(iSR0 < 0)) = 0;


subplot(1,2,1)
plot([0 1.5], [0 1.5], 'k', 'LineWidth', 3);
hold on
plot(extrapA.novel, extrapA.go, '.', 'MarkerSize', 20);
plot(mean(extrapA.novel), mean(extrapA.go), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('M.I. about Go motifs (bits)')
xlabel('M.I. about Novel motifs (bits)')
xlim([0 1.3]);
ylim([0 1.3]);

subplot(1,2,2)
plot([0 3], [0 3], 'k', 'LineWidth', 3);
hold on
plot(extrapAent.novel, extrapAent.go, '.', 'MarkerSize', 20);
plot(mean(extrapAent.novel), mean(extrapAent.go), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('M.I. about Go motifs (bits)')
xlabel('M.I. about Novel motifs (bits)')
xlim([0 3]);
ylim([0 3]);





for i = 1:97,
    text(iSR0(i), iSR1(i), num2str(cellIDs(i)));
end



subplot(1,3,2)
plot([0 1.5], [0 1.5], 'k', 'LineWidth', 3);
hold on
plot(iSR0, iSR2, '.', 'MarkerSize', 20);
plot(mean(iSR0), mean(iSR2), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 14);
ylabel('M.I. about No-Go motifs (bits)')
xlabel('M.I. about Novel motifs (bits)')

for i = 1:97,
    text(iSR0(i), iSR2(i), num2str(cellIDs(i)));
end


subplot(1,3,3)
plot([0 1.5], [0 1.5], 'k', 'LineWidth', 3);
hold on
plot(iSR2, iSR1, '.', 'MarkerSize', 20);
plot(mean(iSR2), mean(iSR1), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 14);
ylabel('M.I. about Go motifs (bits)')
xlabel('M.I. about No-Go motifs (bits)')

for i = 1:97,
    text(iSR2(i), iSR1(i), num2str(cellIDs(i)));
end



%Plot Ent for comparison figure
subplot(2,1,1)
plot([0 3], [0 3], 'k', 'LineWidth', 1);
hold on
plot(extrapAent.novel, extrapAent.go, '.k', 'MarkerSize', 20);
plot(mean(extrapAent.novel), mean(extrapAent.go), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('H_{Go} (bits)')
xlabel('H_{Novel} (bits)')
xlim([0 3]);
ylim([0 3]);

subplot(2,1,2)
plot([0 3], [0 3], 'k', 'LineWidth', 1);
hold on
plot(extrapAent.novel, extrapAent.nogo, '.k', 'MarkerSize', 20);
plot(mean(extrapAent.novel), mean(extrapAent.nogo), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('H_{Nogo} (bits)')
xlabel('H_{Novel} (bits)')
xlim([0 3]);
ylim([0 3]);
set(gcf, 'Position', [0 0 200 400])

%%%%Plot MI for comparison figure
subplot(2,1,1)
plot([0 3], [0 3], 'k', 'LineWidth', 1);
hold on
plot(extrapA.novel, extrapA.go, '.k', 'MarkerSize', 20);
plot(mean(extrapA.novel), mean(extrapA.go), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('MI_{Go} (bits)')
xlabel('MI_{Novel} (bits)')
xlim([0 1.5]);
ylim([0 1.5]);

subplot(2,1,2)
plot([0 3], [0 3], 'k', 'LineWidth', 1);
hold on
plot(extrapA.novel, extrapA.nogo, '.k', 'MarkerSize', 20);
plot(mean(extrapA.novel), mean(extrapA.nogo), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('MI_{Nogo} (bits)')
xlabel('MI_{Novel} (bits)')
xlim([0 1.5]);
ylim([0 1.5]);
set(gcf, 'Position', [0 0 200 400])







subplot(1,2,2)
plot([0 3], [0 3], 'k', 'LineWidth', 3);
hold on
plot(extrapAent.novel, extrapAent.go, '.', 'MarkerSize', 20);
plot(mean(extrapAent.novel), mean(extrapAent.go), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 11);
ylabel('M.I. about Go motifs (bits)')
xlabel('M.I. about Novel motifs (bits)')
xlim([0 3]);
ylim([0 3]);














%Make Information plots highlighting each cell
subplot(1,3,1)
plot([-.1 1.5], [-.1 1.5], 'k', 'LineWidth', 1);
hold on
%plot(iSR0, iSR1, '.', 'MarkerSize', 20);
plot(info_class0(cellIDs ~= 214), info_class1(cellIDs ~= 214), '.k', 'MarkerSize', 20);
plot(info_class0(cellIDs == 214), info_class1(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(info_class0), mean(info_class1), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('MI_{Go} (bits)')
xlabel('MI_{Novel} (bits)')
xlim([-.1 1.5]);
ylim([-.1 1.5]);


subplot(1,3,2)
plot([-.1 1.5], [-.1 1.5], 'k', 'LineWidth', 1);
hold on
%plot(iSR0, iSR2, '.', 'MarkerSize', 20);
plot(info_class0(cellIDs ~= 214), info_class2(cellIDs ~= 214), '.k', 'MarkerSize', 20);
plot(info_class0(cellIDs == 214), info_class2(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(info_class0), mean(info_class2), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('MI_{Nogo} (bits)')
xlabel('MI_{Novel} (bits)')
xlim([-.1 1.5]);
ylim([-.1 1.5]);



subplot(1,3,3)
plot([-.1 1.5], [-.1 1.5], 'k', 'LineWidth', 1);
hold on
plot(info_class2(cellIDs ~= 214), info_class1(cellIDs ~= 214), '.k', 'MarkerSize', 20);
plot(info_class2(cellIDs == 214), info_class1(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(info_class2), mean(info_class1), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('MI_{Go} (bits)')
xlabel('MI_{Nogo} (bits)')
set(gcf, 'Position', [0 0 1000 300]);
xlim([-.1 1.5]);
ylim([-.1 1.5]);




subplot(1,3,1)
plot([0 3], [0 3], 'k', 'LineWidth', 2);
hold on
plot(entR0, entR1, '.', 'MarkerSize', 20);
plot(mean(entR0), mean(entR1), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy of Go motifs (bits)')
xlabel('Entropy of Novel motifs (bits)')

subplot(1,3,2)
plot([0 3], [0 3], 'k', 'LineWidth', 2);
hold on
plot(entR0, entR2, '.', 'MarkerSize', 20);
plot(mean(entR0), mean(entR2), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy of  No-Go motifs (bits)')
xlabel('Entropy about Novel motifs (bits)')

subplot(1,3,3)
plot([0 3], [0 3], 'k', 'LineWidth', 2);
hold on
plot(entR2, entR1, '.', 'MarkerSize', 20);
plot(mean(entR2), mean(entR1), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy about Go motifs (bits)')
xlabel('Entropy about No-Go motifs (bits)')



%Make Plots highlighting example cell

figure
subplot(1,3,1)
plot([0 3], [0 3], 'k', 'LineWidth', 1);
hold on
plot(entR0(cellIDs ~= 214), entR1(cellIDs ~= 214), '.', 'MarkerSize', 20);
plot(entR0(cellIDs == 214), entR1(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(entR0), mean(entR1), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy of Go motifs (bits)')
xlabel('Entropy of Novel motifs (bits)')
set(gcf, 'Position', [0 0 300 300])

subplot(1,3,2)
plot([0 3], [0 3], 'k', 'LineWidth', 2);
hold on
plot(entR0(cellIDs ~= 214), entR2(cellIDs ~= 214), '.', 'MarkerSize', 20);
plot(entR0(cellIDs == 214), entR2(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(entR0), mean(entR2), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy of  No-Go motifs (bits)')
xlabel('Entropy about Novel motifs (bits)')

subplot(1,3,3)
plot([0 3], [0 3], 'k', 'LineWidth', 2);
hold on
plot(entR2(cellIDs ~= 214), entR1(cellIDs ~= 214), '.', 'MarkerSize', 20);
plot(entR2(cellIDs == 214), entR1(cellIDs == 214), '.r', 'MarkerSize', 20);
plot(mean(entR2), mean(entR1), 'gv', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g')
axis square
set(gca, 'FontSize', 14);
ylabel('Entropy about Go motifs (bits)')
xlabel('Entropy about No-Go motifs (bits)')




subplot(1,2,1)
barerror([1;2;3], [mean(class1range); mean(class2range); mean(class0range)], [std(class1range); std(class2range); std(class0range)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('F.R. Range (Hz)');
title('RM ANOVA, p = 0.0525');

class1MinusClass2 = class1range - class2range;
class2MinusClass0 = class2range - class0range;
class1MinusClass0 = class1range - class0range;

subplot(1,2,2)
barerror([1;2;3], [mean(class1MinusClass0); mean(class1MinusClass2); mean(class2MinusClass0)], [std(class1MinusClass0); std(class2MinusClass0); std(class1MinusClass2)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14)
set(gca, 'XTickLabel', {'Go-Nogo', 'Nogo-Novel', 'Go-Novel'});
ylabel('Mean M.I. Difference \pm S.E.M. (bits)')



info_class1(find(info_class1 < 0)) = 0;
info_class2(find(info_class2 < 0)) = 0;
info_class0(find(info_class0 < 0)) = 0;




X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [info_class1(i) 1 i];
    X((i-1)*3+2,:) = [info_class2(i) 2 i];
    X((i-1)*3+3,:) = [info_class0(i) 3 i];
end
RMAOV1(X)
[h,p] = ttest2(info_class1, info_class2)
[h,p] = ttest2(info_class1, info_class0)
[h,p] = ttest2(info_class2, info_class0)



barerror([1;2;3], [mean(info_class1); mean(info_class2); mean(info_class0)], [std(info_class1); std(info_class2); std(info_class0)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Extrapolated Information (Bits)');
title('RM ANOVA, p < 10^{-6}');


bar(1, mean(info_class1), 'g');
hold on
bar(2, mean(info_class2), 'r');
bar(3, mean(info_class0), 'b');
errorbar([1 2 3], [mean(info_class1) mean(info_class2) mean(info_class0)], [std(info_class1) std(info_class2) std(info_class0)]/sqrt(97), '.k', 'MarkerSize', 1);
xlim([0 4]);
set(gca, 'FontSize', 14)
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('MI \pm S.E.M. (bits)');
box off
set(gcf, 'Position', [0 0 300 500]);



class1percent = 100*(info_class1./info_class0)-100;
class2percent = 100*(info_class2./info_class0)-100;



bar(1, mean(class1percent), 'g');
hold on
bar(2, mean(class2percent), 'r');
errorbar([1 2], [mean(class1percent) mean(class2percent)], [std(class1percent) std(class2percent)]/sqrt(97), '.k', 'MarkerSize', 1);
xlim([0 3]);
set(gca, 'FontSize', 14)
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Go', 'Nogo'});
ylabel('Percentage increase in MI over Novel motifs');
box off
set(gcf, 'Position', [0 0 300 400]);




X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [class1max(i) - sponRate(i) 1 i];
    X((i-1)*3+2,:) = [class2max(i) - sponRate(i) 2 i];
    X((i-1)*3+3,:) = [class0max(i) - sponRate(i) 3 i];
end
RMAOV1(X)

X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [sponRate(i) - class1min(i) 1 i];
    X((i-1)*3+2,:) = [sponRate(i) - class2min(i) 2 i];
    X((i-1)*3+3,:) = [sponRate(i) - class0min(i) 3 i];
end
RMAOV1(X)

eRange1 = class1max - sponRate;
eRange2 = class2max - sponRate;
eRange0 = class0max - sponRate;

iRange1 = sponRate - class1min;
iRange2 = sponRate - class2min;
iRange0 = sponRate - class0min;


subplot(1,2,1)
barerror([1;2;3], [mean(eRange1); mean(eRange2); mean(eRange0)], [std(eRange1); std(eRange2); std(eRange0)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Excitatory F.R. Range (Hz)');
title('RM ANOVA, p = 0.034');

subplot(1,2,2)
barerror([1;2;3], -[mean(iRange1); mean(iRange2); mean(iRange0)], [std(iRange1); std(iRange2); std(iRange0)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Inhibitory F.R. Range (Hz)');
title('RM ANOVA, p = 0.072');


offset = 1:2:30;

for i = 1:length(offset),
    normRange1 = class1range - class1min + offset(i);
    normRange2 = class2range - class2min + offset(i);
    normRange0 = class0range - class0min + offset(i);

    dbRange1 = 20*log10(normRange1/offset(i));
    dbRange2 = 20*log10(normRange2/offset(i));
    dbRange0 = 20*log10(normRange0/offset(i));

    meanDBRange1(i) = mean(dbRange1);
    meanDBRange2(i) = mean(dbRange2);
    meanDBRange0(i) = mean(dbRange0);

    stdDBRange1(i) = std(dbRange1);
    stdDBRange2(i) = std(dbRange2);
    stdDBRange0(i) = std(dbRange0);
    
    X = zeros(length(cellIDs)*3,3);
    for j = 1:length(cellIDs),
        X((j-1)*3+1,:) = [dbRange1(j) 1 j];
        X((j-1)*3+2,:) = [dbRange2(j) 2 j];
        X((j-1)*3+3,:) = [dbRange0(j) 3 j];
    end
    p(i) = RMAOV1(X);
    
    
end

subplot(2,1,1)
errorbar(offset,meanDBRange1, stdDBRange1/sqrt(97),'b');
hold on
errorbar(offset,meanDBRange2, stdDBRange2/sqrt(97),'g');
errorbar(offset,meanDBRange0, stdDBRange0/sqrt(97),'r');
xlabel('denominator value');
ylabel('dB of dynamic range');
legend('Go', 'Nogo', 'Novel')

subplot(2,1,2)
plot(offset, p)
xlabel('denominator value')
ylabel('p-value in RM ANOVA')

offset = 1;
normRange1 = class1range - class1min + offset;
normRange2 = class2range - class2min + offset;
normRange0 = class0range - class0min + offset;

dbRange1 = 20*log10(normRange1/offset);
dbRange2 = 20*log10(normRange2/offset);
dbRange0 = 20*log10(normRange0/offset);

subplot(1,2,1)
barerror([1;2;3], [mean(dbRange1); mean(dbRange2); mean(dbRange0)], [std(dbRange1); std(dbRange2); std(dbRange0)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Dynamic Range (dB)');
title('RM ANOVA, p = 0.03');
axis square
ylim([0 30])

subplot(1,2,2)
plot([5 35], [5 35], 'k','LineWidth', 2)
set(gca, 'FontSize', 14);
hold on
plot(dbRange0, dbRange1, '.', 'MarkerSize', 20)
plot(mean(dbRange0), mean(dbRange1), 'ro', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r')

xlim([5 35])
ylim([5 35])
axis square
xlabel('Novel Motif Dynamic Range (dB)')
ylabel('Go Motif Dynamic Range (dB)')



X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [dbRange1(i) 1 i];
    X((i-1)*3+2,:) = [dbRange2(i) 2 i];
    X((i-1)*3+3,:) = [dbRange0(i) 3 i];
end
RMAOV1(X)

[h,p] = ttest(dbRange1-dbRange2)
[h,p] = ttest(dbRange1-dbRange0)
[h,p] = ttest(dbRange2-dbRange0)



X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [class1std(i) 1 i];
    X((i-1)*3+2,:) = [class2std(i) 2 i];
    X((i-1)*3+3,:) = [class0std(i) 3 i];
end
RMAOV1(X)


X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [class1max(i) 1 i];
    X((i-1)*3+2,:) = [class2max(i) 2 i];
    X((i-1)*3+3,:) = [class0max(i) 3 i];
end
RMAOV1(X)

[h,p] = ttest(class1range-class0range)
[h,p] = ttest(class2range-class0range)
[h,p] = ttest(class1range-class2range)


subplot(1,2,1)
barerror([1;2;3], [mean(class1std); mean(class2std); mean(class0std)], [std(class1std); std(class2std); std(class0std)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('std of F.R. distribution (Hz)');
title('RM ANOVA, p = 0.0016');

subplot(1,2,2)
barerror([1;2;3], [mean(class1max); mean(class2max); mean(class0max)], [std(class1max); std(class2max); std(class0max)]/sqrt(97), 1, 'r', 'k')
set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Max F.R. (Hz)');
title('RM ANOVA, p = 0.0341');





subplot(1,2,1)
%barerror([1;2;3], [mean(iSR1); mean(iSR2); mean(iSR0)], [std(iSR1); std(iSR2); std(iSR0)]/sqrt(92), 1, 'r', 'k')
bar(1,mean(iSR1), 'r');
hold on
bar(2,mean(iSR2), 'g');
bar(3,mean(iSR0), 'b');
h = errorbar([1 2 3], mean([iSR1' iSR2' iSR0']), std([iSR1' iSR2' iSR0'])/sqrt(97), '.k');
set(gca, 'FontSize', 14);
set(h, 'MarkerSize', 1);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('MI \pm S.E.M. (bits)');
ylim([0 0.7])
box off
set(gcf, 'Position', [0 0 300 500])

set(gca, 'FontSize', 14);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('I(Motif;FR) (bits)');
title('N = 92 cells, RM ANOVA, p < 0.00002');

class1MinusClass2 = iSR1 - iSR2;
class2MinusClass0 = iSR2 - iSR0;
class1MinusClass0 = iSR1 - iSR0;


pairwise = load('PairwiseInfo');
h = barerror([1;2;3], [mean(iSR1) mean(pairwise.iSR1); mean(iSR2) mean(pairwise.iSR2); mean(iSR0) mean(pairwise.iSR0) ], [[std(iSR1); std(iSR2); std(iSR0)]/sqrt(97) [std(pairwise.iSR1); std(pairwise.iSR2); std(pairwise.iSR0)]/sqrt(597)], 1, 'rb', 'kk');
ylim([0 1.21]);
set(gca, 'FontSize', 14)
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Mean M.I. \pm S.E.M. (bits)')
legend(h, 'Single Unit Info', 'Pairwise Info')


subplot(1,2,2)
barerror([1;2;3], mean([class1MinusClass2' class2MinusClass0' class1MinusClass0'])', std([class1MinusClass2' class2MinusClass0' class1MinusClass0'])'/sqrt(92), 1,'r', 'k')
ylim([0 0.2])
set(gca, 'FontSize', 14)
set(gca, 'XTickLabel', {'Go-Nogo', 'Nogo-Novel', 'Go-Novel'});
ylabel('Mean M.I. Difference \pm S.E.M. (bits)')
title('N = 92 cells')
text(0.9,0.105, 'p < 0.0006')
text(1.85,0.055, 'p < 0.3')
text(2.8,0.14, 'p < 0.00003')


[h,p] = ttest(class1MinusClass2)
[h,p] = ttest(class1MinusClass0)
[h,p] = ttest(class2MinusClass0)



subplot(1,2,1)
barerror([1;2;3], [mean(entR1); mean(entR2); mean(entR0)], [std(entR1); std(entR2); std(entR0)]/sqrt(92), 1, 'r', 'k')
set(gca, 'FontSize', 14)
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Mean Entropy \pm S.E.M. (bits)')

bar(1,mean(entR1), 'r');
hold on
bar(2,mean(entR2), 'g');
bar(3,mean(entR0), 'b');
h = errorbar([1 2 3], mean([entR1' entR2' entR0']), std([entR1' entR2' entR0'])/sqrt(97), '.k');
set(gca, 'FontSize', 14);
set(h, 'MarkerSize', 1);
set(gca, 'XTick', [1 2 3]);
set(gca, 'XTickLabel', {'Go', 'Nogo', 'Novel'});
ylabel('Mean Entropy \pm S.E.M. (bits)')
ylim([0 2.4])
box off
set(gcf, 'Position', [0 0 300 450])



title('N = 92 cells, RM ANOVA p < 0.001')


class1MinusClass2 = entR1 - entR2;
class2MinusClass0 = entR2 - entR0;
class1MinusClass0 = entR1 - entR0;

[h,p] = ttest(entR1-entR2)
[h,p] = ttest(entR1-entR0)
[h,p] = ttest(entR2-entR0)


subplot(1,2,2)
barerror([1;2;3], mean([class1MinusClass2' class2MinusClass0' class1MinusClass0'])', std([class1MinusClass2' class2MinusClass0' class1MinusClass0'])'/sqrt(92), 1,'r', 'k')
ylim([0 0.2])
set(gca, 'FontSize', 14)
set(gca, 'XTickLabel', {'Go-Nogo', 'Nogo-Novel', 'Go-Novel'});
ylabel('Mean Entropy Difference \pm S.E.M. (bits)')
title('N = 92 cells')

text(0.9,0.1, 'p < 0.05')
text(1.8,0.13, 'p < 0.04')
text(2.8,0.17, 'p < 0.001')




X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [entR1(i) 1 i];
    X((i-1)*3+2,:) = [entR2(i) 2 i];
    X((i-1)*3+3,:) = [entR0(i) 3 i];
end
RMAOV1(X)

X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [entR1_raw(i) 1 i];
    X((i-1)*3+2,:) = [entR2_raw(i) 2 i];
    X((i-1)*3+3,:) = [entR0_raw(i) 3 i];
end
RMAOV1(X)




%Run RM anova on the data:

X = zeros(length(cellIDs)*3,3);
for i = 1:length(cellIDs),
    X((i-1)*3+1,:) = [iSR1(i) 1 i];
    X((i-1)*3+2,:) = [iSR2(i) 2 i];
    X((i-1)*3+3,:) = [iSR0(i) 3 i];
end


RMAOV1(X)

plot([1 2 3], [iSR1; iSR2; iSR0])
plot([1 2 3], zscore([iSR1; iSR2; iSR0])')

% ziSR1 = zscore(iSR1);
% ziSR2 = zscore(iSR2);
% ziSR0 = zscore(iSR0);
% 
% 
% for i = 1:length(cellIDs),
%     X((i-1)*3+1,:) = [ziSR1(i) 1 i];
%     X((i-1)*3+2,:) = [ziSR2(i) 2 i];
%     X((i-1)*3+3,:) = [ziSR0(i) 3 i];
% end



barerror([1 2 3]', mean(zscore([iSR1; iSR2; iSR0])')', std(zscore([iSR1; iSR2; iSR0])')'/sqrt(120), 1, 'k', 'k')

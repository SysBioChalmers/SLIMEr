%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-01-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read data and modify for plotting:
addpath('../data')
addpath('../models')
data = readData;

%1. Exp data: lipid classes
lipids     = data.lipidData.metAbbrev([1,3:end]);       %Take out ergosterol
lipidNames = data.lipidData.metNames([1,3:end]);        %Take out ergosterol
abundance  = data.lipidData.abundance([1,3:end])*1000;  %mg/gDW
color      = [0  1  0];                                 %Green
barPlot(abundance,lipids,'[mg/gDW]',color,15,600);

%2. Exp data: chains
chains    = data.chainData.metNames(1:end-2);       %Take out very long chain F.A.s
chains    = strrep(chains,' chain','');
chains    = strrep(chains,'C','');
abundance = data.chainData.abundance(1:end-2)*1000; %mg/gDW
color     = [1  0  0];                              %Red
barPlot(abundance,chains,'[mg/gDW]',color,25,500);

%3. Compare distributions of chains, with same amount of glucose:
new = getLipidDistribution(model_SLIMEr,lipidNames,chains,data.fluxData);
data.fluxData.averages(1) = new.vgluc;
data.fluxData.stdevs(1)   = 0;
old        = getLipidDistribution(model_correctedComp,lipidNames,chains,data.fluxData);
oldChains  = sum(old.comp)';
newChains  = sum(new.comp)';
abundances = [oldChains/sum(oldChains) newChains/sum(newChains)]*100;
color      = [1  1  0      %Yellow
              1  0  0];    %Red
barPlot(abundances,chains,'[%]',color,60,900);
legend('Yeast7 - correct lipid composition','Yeast7 - correct lipid+chain composition','Location','northwest')
legend('boxoff')

%4. Compare distribution of chains in lipids - old model:
color = [0    1    0
         0.2  0.7  0
         0.5  0.5  0
         0.8  0.3  0
         1    0    0];
barPlot(old.comp,lipids,'[mg/gDW]',color,20,900);
legend(chains,'Location','northwest')
legend('boxoff')

%5. Compare variability of chains in lipids - old model:
b = barPlot(old,lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%6. Compare distribution of chains in lipids - new model:
barPlot(new.comp,lipids,'[mg/gDW]',color,20,900);
legend(chains,'Location','northwest')
legend('boxoff')

%7. Compare variability of chains in lipids - new model:
b = barPlot(new,lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Compare energy differences:
netATP  = num2str(round(old.netATP - new.netATP,2));
GAMunk  = 35.01;    %mmol/gDW (old model from createNewModels.m)
percATP = num2str(round((old.netATP - new.netATP)/GAMunk*100,1));
disp(['Net ATP spent in changing lipid comp: ' netATP ' mmol/gDW = ' percATP '% of GAM'])
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-01-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Gety lipid distribution for each condition:
addpath('../data')
addpath('../models')
data = cell(size(model_correctedComp));
old  = cell(size(model_correctedComp));
new  = cell(size(model_correctedComp));
for i = 1:length(data)
    %Read and modify data:
    data{i}    = readData(i);
    lipidNames = data{i}.lipidData.metNames([1,3:end]);     %Take out ergosterol
    chains     = data{i}.chainData.metNames(1:end-2);       %Take out very long chain F.A.s
    chains     = strrep(chains,' chain','');
    chains     = strrep(chains,'C','');
    
    % Compute distributions of chains, with same amount of glucose:
    new{i} = getLipidDistribution(model_SLIMEr{i},lipidNames,chains,data{i}.fluxData);
    data{i}.fluxData.averages(1) = new{i}.vgluc;
    data{i}.fluxData.stdevs(1)   = 0;
    old{i} = getLipidDistribution(model_correctedComp{i},lipidNames,chains,data{i}.fluxData);
end

%1. Exp data: lipid classes
lipids     = data{1}.lipidData.metAbbrev([1,3:end]);        %Take out ergosterol
abundance  = data{1}.lipidData.abundance([1,3:end])*1000;	%mg/gDW
color      = [0  1  0];                                     %Green
barPlot(abundance,lipids,'[mg/gDW]',color,15,600);

%2. Exp data: chains
abundance = data{1}.chainData.abundance(1:end-2)*1000;  %mg/gDW
color     = [1  0  0];                                  %Red
barPlot(abundance,chains,'[mg/gDW]',color,25,500);

%3. Compare distributions of chains:
oldChains  = sum(old{1}.comp)';
newChains  = sum(new{1}.comp)';
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
barPlot(old{1}.comp,lipids,'[mg/gDW]',color,20,900);
legend(chains,'Location','northwest')
legend('boxoff')

%5. Compare variability of chains in lipids - old model:
b = barPlot(old{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%6. Compare distribution of chains in lipids - new model:
barPlot(new{1}.comp,lipids,'[mg/gDW]',color,20,900);
legend(chains,'Location','northwest')
legend('boxoff')

%7. Compare variability of chains in lipids - new model:
b = barPlot(new{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Compare energy differences at reference condition:
netATP  = num2str(round(old{1}.netATP - new{1}.netATP,2));
NGAM    = 0.7;
GAMunk  = 42.12 - NGAM/0.1;    %mmol/gDW (old model from createNewModels.m)
percATP = num2str(round((old{1}.netATP - new{1}.netATP)/GAMunk*100,1));
disp(['Net ATP spent in changing lipid comp: ' netATP ' mmol/gDW = ' percATP '% of GAM'])
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
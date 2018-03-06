%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-03-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get lipid distribution for each condition:
close all
addpath('../data')
addpath('../models')
data = cell(size(model_correctedComp));
old  = cell(size(model_correctedComp));
new  = cell(size(model_correctedComp));
for i = 1:length(data)
    %Read and modify data:
    data{i}    = readLahtveeData(i);
    lipidNames = data{i}.lipidData.metNames([1,3:end]);     %Take out ergosterol
    chains     = data{i}.chainData.metNames(1:end-2);       %Take out very long chain F.A.s
    chains     = strrep(chains,' chain','');
    chains     = strrep(chains,'C','');
    
    %Compute distributions of chains, with same amount of glucose:
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

%8. Compare differences at stress conditions:
Cchange = zeros(length(model_SLIMEr),1);
NGAMs   = zeros(length(model_SLIMEr),1);
Cdiff   = zeros(length(model_SLIMEr),1);    
ATPdiff = zeros(length(model_SLIMEr),1);
meanVar = zeros(length(model_SLIMEr),1);    
for i = 1:length(model_SLIMEr)
    %Difference in chain composition to REF [%]:
    refRelComp = sum(new{1}.comp)/sum(sum(new{1}.comp));
    newRelComp = sum(new{i}.comp)/sum(sum(new{i}.comp));
    Cchange(i) = mean(abs(newRelComp - refRelComp))*100;	%%
    
    %Fitted maintenance [mmol/gDWh]:
    NGAMref  = 0.7;      %mmol/gDWh (Nilsson et al. 2016)
    GAMunk   = old{1}.netATP - GAMpol(1) - NGAMref/old{1}.mu;	%mmol/gDW
    NGAMs(i) = (old{i}.netATP - GAMpol(i) - GAMunk)*old{i}.mu;	%mmol/gDWh
    
    %Extra carbon cost [mmol/gDW]:
    formulas = data{i}.chainData.formulas(1:4);
    MWs      = getMWfromFormula(formulas)*1000;                     %g/mol
    Cmol     = cellfun(@(x) x(2:3),formulas,'UniformOutput',false);	%C-mol/mol
    Cmol     = str2double(Cmol);                                    %C-mol/mol
    massDiff = sum(new{i}.comp - old{i}.comp)'/1000;              	%g/gDW
    CmolDiff = massDiff./MWs.*Cmol;                                 %C-mol/gDW
    Cdiff(i) = sum(abs(CmolDiff))*1000;                             %mmol/gDW
    
    %Extra ATP cost [mmol/gDW]
    ATPdiff(i) = old{i}.netATP - new{i}.netATP;	%mmol/gDW
    
    %Mean variability [mg/gDW]
    newVar     = new{i}.var.max - new{i}.var.min;   %mg/gDW
    meanVar(i) = mean(mean(newVar));                %mg/gDW
end
stressData = [Cchange NGAMs meanVar Cdiff ATPdiff];
varNames = {'Experimental difference in chain composition compared to reference [%]', ...
            'Fitted maintenance in original model [mmol/gDWh]', ...
            'Mean variability in model with SLIME rxns [mg/gDW]', ...
            'Extra carbon cost of adding SLIME rxns [mmol/gDW]', ...
            'Extra ATP cost of adding SLIME rxns [mmol/gDW]'};
figure('position', [100,100,900,400])
subplot(1,3,1)
stressPlot([30 33 36 38],stressData(1:4,:),{'','','','',''},'Temperature [°C]')
subplot(1,3,2)
stressPlot([0 200 400 600],stressData([1 5:7],:),varNames,'NaCl [mM]')
subplot(1,3,3)
stressPlot([0 20 40 60],stressData([1 8:10],:),{'','','','',''},'Ethanol [g/L]')
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
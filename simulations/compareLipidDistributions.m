%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-05-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get lipid distribution for each condition:
close all
addpath('../data')
addpath('../models')
data = cell(size(model_corrComp));
old  = cell(size(model_corrComp));
new  = cell(size(model_corrComp));
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
    old{i} = getLipidDistribution(model_corrComp{i},lipidNames,chains,data{i}.fluxData);
end

%Fig 2A: Compare distributions of chains
expChains = data{1}.chainData.abundance(1:end-2)*1000;  %mg/gDW
stdChains = data{1}.chainData.std(1:end-2)*1000;        %mg/gDW
oldChains = sum(old{1}.comp)';
newChains = sum(new{1}.comp)';
oldChains = oldChains/sum(oldChains)*100;
newChains = newChains/sum(newChains)*100;
stdChains = stdChains/sum(expChains)*100;
expChains = expChains/sum(expChains)*100;
color     = [1  1  0      %Yellow
             0  0  1      %Blue
             1  0  0];    %Red
b = barPlot([oldChains newChains expChains],chains,'[%]',color,60,900);
hold on
posExp = b(3).XData + b(3).XOffset;
e = errorbar(posExp',expChains,stdChains,'k.');
e.Marker = 'none';
legend(b(1:3),'Permissive model','Enhanced model', ...
              'Experimental values','Location','northeast')
legend('boxoff')
hold off

%Fig S2: Compare variability of chains in lipids - old model:
lipids = data{1}.lipidData.metAbbrev([1,3:end]);        %Take out ergosterol
color  = [0    1    0
          0.2  0.7  0
          0.5  0.5  0
          0.8  0.3  0
          1    0    0];
b = barPlot(old{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Fig 2B: Compare variability of chains in lipids - new model:
b = barPlot(new{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Compare differences at stress conditions:
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
    ATPdiff(i) = (old{i}.netATP - new{i}.netATP)*1000;  %umol/gDW
    
    %Mean variability [mg/gDW]
    newVar     = new{i}.var.max - new{i}.var.min;       %mg/gDW
    meanVar(i) = mean(mean(newVar));                    %mg/gDW
end

%Fig 3C: Main variables
mainStressData = [Cdiff ATPdiff];
mainVarNames   = {'Extra carbon cost [mmol/gDW]', ...
                  'Extra ATP cost [\mumol/gDW]'};
stressPlot(mainStressData,mainVarNames,[1 3],[50 300])

%Fig S4: Supplementary variables
suppStressData = [Cchange NGAMs meanVar];
suppVarNames   = {'Experimental difference in chain composition compared to reference [%]', ...
                  'Fitted maintenance in permissive model [mmol/gDWh]', ...
                  'Mean variability in enhanced model [mg/gDW]'};
stressPlot(suppStressData,suppVarNames,[0 20],[3 5])
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
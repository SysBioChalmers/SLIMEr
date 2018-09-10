%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamin J. Sanchez. Last update: 2018-09-04
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
    fluxData   = data{i}.fluxData;
    
    %Compute distributions of chains, with same amount of glucose:
    if i == 1
        getFullVar = true;
    else
        getFullVar = false;
    end
    new{i} = getLipidDistribution(model_SLIMEr{i},lipidNames,chains,fluxData,getFullVar);
    data{i}.fluxData.averages(1) = new{i}.vgluc;
    data{i}.fluxData.stdevs(1)   = 0;
    old{i} = getLipidDistribution(model_corrComp{i},lipidNames,chains,fluxData,getFullVar);
    disp(['Computing lipid distribution: ready with dataset ' num2str(i)])
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

%Fig S3: Compare variability of chains in lipids - old model:
lipids = data{1}.lipidData.metAbbrev([1,3:end]);        %Take out ergosterol
color  = sampleCVDmap(4);
b = barPlot(old{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Fig 2B: Compare variability of chains in lipids - new model:
b = barPlot(new{1},lipids,'[mg/gDW]',color,20,900);
legend(b,chains,'Location','northwest')
legend('boxoff')

%Fig S4A: Compare full FVA for reference conditions:
FVA = zeros(length(model_SLIMEr{1}.rxns),2);
for i = 1:length(model_SLIMEr{1}.rxns)
    cross_pos = strcmp(model_corrComp{1}.rxnNames,model_SLIMEr{1}.rxnNames{i});
    if sum(cross_pos) > 0 %take out pseudorxns or SLIMEr stuff
        if sum(cross_pos) > 1
            cross_pos = i;
        end
        FVA(i,1) = old{1}.fullVar(cross_pos);
        FVA(i,2) = new{1}.fullVar(i);
    end
end
FVA = FVA(sum(FVA,2) > 0,:);
figure('position', [100,100,500,500])
hold on
plot([1e-4,1e4],[1e-4,1e4],'-k')
plot(FVA(:,1),FVA(:,2),'ob')
plotOptions([1e-4,1e4],[1e-4,1e4],'Permissive model variability [mmol/gDWh]', ...
    'Enhanced model variability [mmol/gDWh]','','','','',12)
set(gca,'YScale','log')
set(gca,'XScale','log')
hold off

%Compare differences at stress conditions:
Cchange = zeros(length(model_SLIMEr),1);
NGAMs   = zeros(length(model_SLIMEr),1);
Cdiff   = zeros(length(model_SLIMEr),1);    
ATPdiff = zeros(length(model_SLIMEr),1);
lipVar  = zeros(length(model_SLIMEr),1);    
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
    newVar    = new{i}.lipVar.max - new{i}.lipVar.min;	%mg/gDW
    lipVar(i) = mean(mean(newVar));                     %mg/gDW
end

%Fig 3C: Main variables
mainStressData = [Cdiff ATPdiff];
mainVarNames   = {'Extra carbon cost [mmol/gDW]', ...
                  'Extra ATP cost [\mumol/gDW]'};
stressPlot(mainStressData,mainVarNames,[1 3],[50 300])

%Fig S6: Supplementary variables
suppStressData = [Cchange NGAMs lipVar];
suppVarNames   = {'Experimental difference in chain composition compared to reference [%]', ...
                  'Fitted maintenance in permissive model [mmol/gDWh]', ...
                  'Mean variability in enhanced model [mg/gDW]'};
stressPlot(suppStressData,suppVarNames,[0 20],[3 5])
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

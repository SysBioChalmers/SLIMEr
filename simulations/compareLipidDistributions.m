%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-01-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. Exp data: lipid classes
fid = fopen('../data/lipid_data.csv');
lipidData  = textscan(fid,'%s %s %s %f32','Delimiter',',','HeaderLines',1);
lipids     = lipidData{1}([1,3:end]);       %Take out ergosterol
lipidNames = lipidData{2}([1,3:end]);       %Take out ergosterol
data       = lipidData{4}([1,3:end])*1000;  %mg/gDW
fclose(fid);
color = [0  1  0];    %Green
barPlot(data,lipids,'[mg/gDW]',color,15,600);

%2. Exp data: chains
fid = fopen('../data/chain_data.csv');
chainData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
chains    = chainData{1}(1:end-2);      %Take out very long chain F.A.s
chains    = strrep(chains,' chain','');
chains    = strrep(chains,'C','');
data      = chainData{3}(1:end-2)*1000;	%mg/gDW
fclose(fid);
color = [1  0  0];    %Red
barPlot(data,chains,'[mg/gDW]',color,25,500);

%3. Compare distributions of chains:
old       = getLipidDistribution(model_correctedComp,lipidNames,chains);
new       = getLipidDistribution(model_SLIMEr,lipidNames,chains);
oldChains = sum(old.comp)';
newChains = sum(new.comp)';
data      = [oldChains/sum(oldChains) newChains/sum(newChains)]*100;
color     = [1  1  0      %Yellow
             1  0  0];    %Red
barPlot(data,chains,'[%]',color,60,900);
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
netATP  = num2str(round(new.netATP - old.netATP,2));
ATPpos  = strcmp(model_correctedComp.mets,'s_0434[c]');
Xpos    = strcmp(model_correctedComp.rxns,'r_4041');
GAM     = abs(model_correctedComp.S(ATPpos,Xpos));
percATP = num2str(round((new.netATP - old.netATP)/GAM*100,1));
disp(['net ATP spent in changing lipid comp: ' netATP ' mmol/gDW = ' percATP '% of GAM'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
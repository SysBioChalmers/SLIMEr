%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2017-12-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1. Exp data: lipid classes
cd ../data
fid = fopen('lipid_data.csv');
lipidData  = textscan(fid,'%s %s %s %f32','Delimiter',',','HeaderLines',1);
lipids     = lipidData{1}([1,3:end]);       %Take out ergosterol
lipidNames = lipidData{2}([1,3:end]);       %Take out ergosterol
data       = lipidData{4}([1,3:end])*1000;  %mg/gDW
fclose(fid);
cd ../simulations
color = [0  1  0];    %Green
barPlot(data,lipids,'[mg/gDW]',color,15,600)

%2. Exp data: chains
cd ../data
fid = fopen('chain_data.csv');
chainData = textscan(fid,'%s %s %f32','Delimiter',',','HeaderLines',1);
chains    = chainData{1}(1:end-2);      %Take out very long chain F.A.s
chains    = strrep(chains,' chain','');
chains    = strrep(chains,'C','');
data      = chainData{3}(1:end-2)*1000;	%mg/gDW
fclose(fid);
cd ../simulations
color = [1  0  0];    %Red
barPlot(data,chains,'[mg/gDW]',color,25,500)

%3. Compare distributions of chains:
[comp_old,~]   = getLipidDistribution(model_correctedComp,lipidNames,chains);
comp_old       = comp_old*1000;     %mg/gDW
[comp_new,var] = getLipidDistribution(model_SLIMEr,lipidNames,chains);
comp_new       = comp_new*1000;     %mg/gDW
chains_old     = sum(comp_old)';
chains_new     = sum(comp_new)';
data           = [chains_old/sum(chains_old) chains_new/sum(chains_new)]*100;
color          = [1  1  0      %Yellow
                  1  0  0];    %Red
barPlot(data,chains,'[%]',color,60,900)
legend('Yeast7 - correct lipid composition','Yeast7 - correct lipid+chain composition','Location','northwest')
legend('boxoff')

%4. Compare distribution of chains in lipids:
color = [0    1    0
         0.2  0.7  0
         0.5  0.5  0
         0.8  0.3  0
         1    0    0];
barPlot(comp_new,lipids,'[mg/gDW]',color,10,900)
legend(chains,'Location','northwest')
legend('boxoff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
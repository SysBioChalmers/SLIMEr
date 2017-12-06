%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions(model_correctedComp,model_SLIMEr)
%
% Benjamín J. Sánchez. Last update: 2017-12-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compareLipidDistributions(model_correctedComp,model_SLIMEr)

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
data_old   = getLipidDistribution(model_correctedComp,lipidNames,chains);
data_old   = data_old*1000;     %mg/gDW
data_new   = getLipidDistribution(model_SLIMEr,lipidNames,chains);
data_new   = data_new*1000;     %mg/gDW
chains_old = sum(data_old)';
chains_new = sum(data_new)';
data       = [chains_old/sum(chains_old) chains_new/sum(chains_new)]*100;
color      = [1  1  0      %Yellow
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
barPlot(data_new,lipids,'[mg/gDW]',color,5,900)
legend(chains,'Location','northwest')
legend('boxoff')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function barPlot(data,names,units,color,ymax,xlength)

%Plot data:
figure('position', [100,100,xlength,400])
hold on
b = bar(data,'BarWidth',1);

for i = 1:length(b)
    set(b(i),'FaceColor',color(i,:));
end

%Various options:
text_size = 12;
set(gca,'XTick',1:length(names),'XTickLabel',names)
set(gca,'FontSize',text_size)
ylabel(['Abundance ' units],'FontSize',text_size);
ylim([0 ymax])
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
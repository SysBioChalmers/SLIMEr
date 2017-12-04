%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2017-12-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%chains:
chains = {'C14:0','C16:0','C16:1','C18:0','C18:1'};

%Compare distributions of chains:
chains_old = getChainDistribution(model_correctedComp,chains);
chains_new = getChainDistribution(model_SLIMEr,chains);
data       = [chains_old/sum(chains_old) chains_new/sum(chains_new)]*100;

%Plot data:
figure('position', [100,100,900,400])
hold on
b = bar(data,'BarWidth',1);
color = [1  1  0      %Yellow
         1  0  0];    %Red
legend(b,'Yeast7 - correct lipid composition','Yeast7 - correct lipid+chain composition','Location','northwest')

for i = 1:length(b)
    set(b(i),'FaceColor',color(i,:));
end

%Various options:
text_size = 12;
set(gca,'XTick',1:5,'XTickLabel',chains)
set(gca,'FontSize',text_size)
ylabel('Abundance [%]','FontSize',text_size);
legend('boxoff')
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
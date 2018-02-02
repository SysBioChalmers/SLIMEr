%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stressPlot(stressLevels,stressData,varNames,x_lab)
%
% Benjamín J. Sánchez. Last update: 2018-02-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stressPlot(stressLevels,stressData,varNames,x_lab)

cmap = colormap;
step = floor(length(cmap)/4);
hold on
yyaxis left
plot(stressLevels,stressData(:,1),'-ob','Color',cmap(1,:),'LineWidth',2)
plot(stressLevels,stressData(:,2),'-or','Color',cmap(step,:),'LineWidth',2)
plot(stressLevels,stressData(:,3),'-og','Color',cmap(2*step,:),'LineWidth',2)
yyaxis right
plot(stressLevels,stressData(:,4),'-oy','Color',cmap(3*step,:),'LineWidth',2)
plot(stressLevels,stressData(:,5),'-om','Color',cmap(4*step,:),'LineWidth',2)
legend(varNames,'Location','southoutside')
legend('boxoff')
if isempty(varNames{1})
    legend('hide')
end

%Various options:
text_size = 12;
set(gca,'XTick',stressLevels)
xlabel(x_lab,'FontSize',text_size);
xlim([stressLevels(1) stressLevels(end)])
ylim([0 3])
set(gca,'YTick',0:3)
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
yyaxis left
ylim([0 20])
set(gca,'YTick',0:5:20)
set(gca,'FontSize',text_size)
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
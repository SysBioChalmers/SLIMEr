%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stressPlot(stressData,varNames,ylim1,ylim2)
%
% Benjamín J. Sánchez. Last update: 2018-04-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stressPlot(stressData,varNames,ylim1,ylim2)
[~,m] = size(stressData);
if m == 2
    figure('position', [100,100,500,600])
    varNames = {varNames{1},varNames{1},varNames{1},varNames{2},varNames{2},varNames{2}};
    stressSubPlot(0:3,stressData(1:4,:),varNames,'Temperature [°C]',ylim1,ylim2)
    stressSubPlot(0:3,stressData([1 5:7],:),varNames,'NaCl [mM]',ylim1,ylim2)
    stressSubPlot(0:3,stressData([1 8:10],:),varNames,'Ethanol [g/L]',ylim1,ylim2)
else
    figure('position', [100,100,900,500])
    subplot(1,3,1)
    stressSubPlot([30 33 36 38],stressData(1:4,:),{'','',''},'Temperature [°C]',ylim1,ylim2)
    subplot(1,3,2)
    stressSubPlot([0 200 400 600],stressData([1 5:7],:),varNames,'NaCl [mM]',ylim1,ylim2)
    subplot(1,3,3)
    stressSubPlot([0 20 40 60],stressData([1 8:10],:),{'','',''},'Ethanol [g/L]',ylim1,ylim2)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stressSubPlot(stressLevels,stressData,varNames,x_lab,ylim1,ylim2)

if contains(x_lab,'Temp')
    color = [1 0 0];            %red
elseif contains(x_lab,'NaCl')
    color = [0 0 1];            %blue
elseif contains(x_lab,'Ethanol')
    color = [0 128/255 0];      %forest green
end
hold on
yyaxis left
[~,m] = size(stressData);
plot(stressLevels,stressData(:,1),'-o','Color',color,'LineWidth',2)
if m == 3
    plot(stressLevels,stressData(:,2),':o','Color',color,'LineWidth',2)
end
yyaxis right
plot(stressLevels,stressData(:,m),'--o','Color',color,'LineWidth',2)
legend(varNames,'Location','southoutside')
legend('boxoff')
if isempty(varNames{1})
    legend('hide')
end

%Various options:
text_size = 12;
set(gca,'XTick',stressLevels)
if m == 2
    xlabel('Stress level','FontSize',text_size);
else
    xlabel(x_lab,'FontSize',text_size);
end
xlim([stressLevels(1) stressLevels(end)])
ylim(ylim2)
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
yyaxis left
ylim(ylim1)
set(gca,'FontSize',text_size)
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
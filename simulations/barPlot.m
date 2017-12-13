%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = barPlot(data,names,units,color,ymax,xlength)
%
% Benjamín J. Sánchez. Last update: 2017-12-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = barPlot(data,names,units,color,ymax,xlength)

figure('position', [100,100,xlength,400])
hold on

if isfield(data,'var')
    %Variability plot:
    [M,N] = size(data.comp);
    for i = 1:M
        for j = 1:N
            comp  = data.comp(i,j);
            min   = data.var.min(i,j);
            max   = data.var.max(i,j);
            left  = i - 0.5 + (j-1)/N;
            right = i - 0.5 + j/N;
            b(j) = patch([left,right,right,left],[min,min,max,max],color(j,:));
            plot([left,right],[comp,comp],'-k','LineWidth',1.5)
        end
        plot([right,right],[0,ymax],'--k')
    end
else
    %Simple barplot:
    [M,N] = size(data);
    b = bar(data,'BarWidth',1);
    for i = 1:N
        set(b(i),'FaceColor',color(i,:));
    end
end

%Various options:
text_size = 12;
set(gca,'XTick',1:length(names),'XTickLabel',names)
set(gca,'FontSize',text_size)
ylabel(['Abundance ' units],'FontSize',text_size);
xlim([0.5,M+0.5])
ylim([0 ymax])
set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
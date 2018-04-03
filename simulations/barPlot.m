%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = barPlot(data,names,units,color,ymax,xlength,std)
%
% Benjamín J. Sánchez. Last update: 2018-04-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = barPlot(data,names,units,color,ymax,xlength,std)

if nargin < 7
    figure('position', [100,100,xlength,400])
    std = [];
end
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

%Add standard deviation:
if ~isempty(std)
    errorbar(1:length(names),data,std,'.')
end

%Various options:
text_size = 12;
set(gca,'XTick',1:length(names),'XTickLabel',names)
set(gca,'FontSize',text_size)
ylabel(['Abundance ' units],'FontSize',text_size);
xlim([0.5,M+0.5])
ylim([0 ymax])
set(gca,'XColor','k')
set(gca,'YColor','k')
box on
hold off

%Shrink and rotate labels when there are too many:
if length(names) > 40
    axis = gca;
    set(axis,'FontSize',text_size/2)
    xtickangle(90)
    axis.XAxis.TickLength = [0,0];
    if xlength ~= 1400
        set(axis,'XTick',[])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = barPlot(data,names,units,color,ymax,xlength,std,stacked)
%
% Benjamin J. Sanchez. Last update: 2018-09-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = barPlot(data,names,units,color,ymax,xlength,std,stacked)

if nargin < 7
    figure('position', [100,100,xlength,400])
    std = [];
end
if nargin < 8
    stacked = false;
end
hold on

if isfield(data,'lipVar')
    %Variability plot:
    [M,N] = size(data.comp);
    for i = 1:M
        for j = 1:N
            comp  = data.comp(i,j);
            min   = data.lipVar.min(i,j);
            max   = data.lipVar.max(i,j);
            left  = i - 0.5 + (j-1)/N;
            right = i - 0.5 + j/N;
            b(j) = patch([left,right,right,left],[min,min,max,max],color(j,:));
            plot([left,right],[comp,comp],'-k','LineWidth',2.5)
        end
        plot([right,right],[0,ymax],'--k')
    end
else
    %Simple barplot:
    [M,N] = size(data);
    if stacked
        b = bar(data,'stacked','BarWidth',0.6);
    else
        b = bar(data,'BarWidth',1);
    end
    for i = 1:N
        set(b(i),'FaceColor',color(i,:));
    end
end

%Add standard deviation:
if ~isempty(std)
    e = errorbar(1:length(names),data,std,'k.');
    e.Marker = 'none';
end

%Various options:
text_size = 12;
plotOptions([0.5,M+0.5],[0 ymax],'',['Abundance ' units],1:length(names),'',names,'',text_size)
axis = gca;
axis.XAxis.TickLength = [0,0];
hold off

%Shrink and rotate labels when there are too many:
if length(names) > 40
    set(axis,'FontSize',text_size/2)
    xtickangle(90)
    if xlength ~= 1400
        set(axis,'XTick',[])
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

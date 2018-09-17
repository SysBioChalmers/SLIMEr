%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotOptions(x_lim,y_lim,x_lab,y_lab,x_tick,y_tick,x_ticklab,y_ticklab,text_size)
%
% Benjamin J. Sanchez. Last update: 2018-09-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotOptions(x_lim,y_lim,x_lab,y_lab,x_tick,y_tick,x_ticklab,y_ticklab,text_size)

axis = gca;

if ~isempty(x_lim)
    xlim(x_lim)
end

if ~isempty(y_lim)
    ylim(y_lim)
end

if ~isempty(x_lab)
    xlabel(x_lab,'FontSize',text_size);
end

if ~isempty(y_lab)
    ylabel(y_lab,'FontSize',text_size);
end

if ~isempty(x_tick)
    axis.XAxis.TickValues = x_tick;
end

if ~isempty(y_tick)
    axis.YAxis.TickValues = y_tick;
end

if ~isempty(x_ticklab)
    axis.XAxis.TickLabel = x_ticklab;
end

if ~isempty(y_ticklab)
    axis.YAxis.TickLabel = y_ticklab;
end

set(gca,'FontSize',text_size)
set(gca,'XColor','k')
set(gca,'YColor','k')
box on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-04-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
close all
addpath('../data')
addpath('../models')

conds  = {'BY4741 wild type - 24ºC','BY4741 wild type - 37ºC', ...
          'BY4741 \itelo1\rm\Delta - 24ºC','BY4741 \itelo1\rm\Delta - 37ºC', ...
          'BY4741 \itelo2\rm\Delta - 24ºC','BY4741 \itelo2\rm\Delta - 37ºC', ...
          'BY4741 \itelo3\rm\Delta - 24ºC','BY4741 \itelo3\rm\Delta - 37ºC'};
errors = zeros(2,length(model_corrComp_val));
for i = 1:length(model_corrComp_val)
    %Read experimental data for each condition:
    data = readEjsingData(i);
    data = convertEjsingData(data,model,false);
    data.metNames  = data.metNames(1:end-1);         %filter out ergosterol
    data.abundance = data.abundance(1:end-1)*1000;   %mg/gDW
    data.std       = data.std(1:end-1)*1000;         %mg/gDW
    
    Nsim = 10000;
    abundance_modC = lipidRandomSampling(model_corrComp_val{i},data,Nsim);
    abundance_modS = lipidRandomSampling(model_SLIMEr_val{i},data,Nsim);
    
    %Averages & errors:
    means       = [mean(abundance_modC,2) mean(abundance_modS,2)];
    diffs       = means - ones(size(means)).*data.abundance;
    errors(:,i) = sum(abs(diffs'),2)/length(diffs);
    
    %Reshape model results:
    x = repmat((1:length(data.metNames))',Nsim,1);
    y = reshape(abundance_modS,[length(data.metNames)*Nsim,1]);
    
    %Plot data:
    ymax = ceil(max(y)/5)*5;
    if i == 1
        %Fig 3: Random sampling at reference conditions
        xlength = 1400;
        figure('position', [100,100,xlength,400])
        barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
        hold on
        scatter(x,y,10,'b','MarkerEdgeAlpha',.02)
        hold off
        figure('position', [100,100,xlength,600])
    else
        xlength = 1000;
    end
    %Fig S3: Random sampling at all conditions
    subplot(length(model_corrComp_val),1,length(model_corrComp_val)+1-i)
    barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
    hold on
    scatter(x,y,10,'b','MarkerEdgeAlpha',.02)
    text(80,ymax-5,conds{i})
    hold off
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
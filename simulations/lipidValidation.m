%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-04-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
close all
addpath('../data')
addpath('../models')

errors = zeros(length(model_corrComp_val),2);
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
    errors(i,:) = sum(abs(diffs),1)/length(diffs);
    
    %Reshape model results:
    x = repmat((1:length(data.metNames))',Nsim,1);
    y = reshape(abundance_modS,[length(data.metNames)*Nsim,1]);
    
    %Plot data:
    ymax = ceil(max(y)/5)*5;
    barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,1400,data.std);
    hold on
    scatter(x,y,10,'b','MarkerEdgeAlpha',.02)
    hold off
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
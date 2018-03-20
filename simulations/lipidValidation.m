%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get experimental data for each condition:
close all
addpath('../data')
addpath('../models')
data = cell(1,8);
for i = 1:length(data)
    %Read data:
    data = readEjsingData(i);
    data = convertEjsingData(data,model);
    
    %Filter out ergosterol:
    metNames  = data.metNames(1:end-1);
    abundance = data.abundance(1:end-1)*1000;   %mg/gDW
    std       = data.std(1:end-1)*1000;         %mg/gDW
    
    %Plot data:
    barPlot(abundance,metNames,'[mg/gDW]','r',25,1400,std);
    
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
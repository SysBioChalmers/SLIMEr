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
    data = readEjsingData(i,model);
        
    %Filter out ergosterol:
    metNames  = data.metNames(1:end-1);
    abundance = data.abundance(1:end-1);
    std       = data.std(1:end-1);
    
    %Plot data:
    barPlot(abundance,metNames,'[% g/gDW]','r',20,1400,std)
    
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
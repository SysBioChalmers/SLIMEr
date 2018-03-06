%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-03-06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get experimental data for each condition:
close all
addpath('../data')
addpath('../models')
data = cell(1,8);
for i = 1:length(data)
    %Read data:
    data = readEjsingData(i);
        
    %Filter out ergosterol:
    metNames  = data.metNames(1:end-1);
    abundance = data.abundance(1:end-1);
    std       = data.std(1:end-1);
    
    %Plot data:
    barPlot(abundance,metNames,'[% mol/mol]','r',20,1400,std)
    
end
rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data,scaling)
%
% Benjamín J. Sánchez. Last update: 2018-03-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data,scaling)

%Model with lipid composition corrected:
model_corrComp = SLIMEr(model,data,false);

%Model with both lipid and chain length constrained to data:
model_SLIMEr = SLIMEr(model,data,true);

%Make abundances be consistent:
[model_SLIMEr,k] = scaleAbundancesInModel(model_SLIMEr,data,scaling);
model_corrComp   = adjustModel(model_corrComp,k,false,scaling);

%Correct the rest of the composition to be consistent:
[model_corrComp,~]    = changeOtherComp(model_corrComp,data);
[model_SLIMEr,GAMpol] = changeOtherComp(model_SLIMEr,data);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
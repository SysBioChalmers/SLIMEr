%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data)
%
% Benjamín J. Sánchez. Last update: 2018-03-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data)

%Model with lipid composition corrected:
model_corrComp = SLIMEr(model,data,false);

%Model with both lipid and chain length constrained to data:
model_SLIMEr = SLIMEr(model,data,true);

%Make abundances be consistent:
[model_SLIMEr,k] = scaleAbundancesInModel(model_SLIMEr,data);
model_corrComp   = adjustModel(model_corrComp,k,false);

%Correct the rest of the composition to be consistent:
[model_corrComp,~]    = changeOtherComp(model_corrComp,data);
[model_SLIMEr,GAMpol] = changeOtherComp(model_SLIMEr,data);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
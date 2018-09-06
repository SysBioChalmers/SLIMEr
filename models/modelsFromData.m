%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data,scaling)
%
% Benjamin J. Sanchez. Last update: 2018-09-05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model_corrComp,model_SLIMEr,k,GAMpol] = modelsFromData(model,data,scaling)

%Model with lipid composition corrected:
model_corrComp = SLIMEr(model,data,false);

%Model with both lipid and chain length constrained to data:
model_SLIMEr = SLIMEr(model,data,true);

%Correct the biomass composition with data:
[model_corrComp,~] = changeOtherComp(model_corrComp,data);
[model_SLIMEr,~]   = changeOtherComp(model_SLIMEr,data);

%Make abundances of backbones & chains be consistent:
[model_SLIMEr,k] = scaleAbundancesInModel(model_SLIMEr,data,scaling);
model_corrComp   = adjustModel(model_corrComp,k,false,scaling);

%Correct again the biomass composition, to fix the new lipid content:
[model_corrComp,~]    = changeOtherComp(model_corrComp,data);
[model_SLIMEr,GAMpol] = changeOtherComp(model_SLIMEr,data);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

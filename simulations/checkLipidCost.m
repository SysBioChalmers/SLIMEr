%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compareLipidDistributions
%
% Benjamín J. Sánchez. Last update: 2018-04-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function costs = checkLipidCost(model)

addpath('../data')
data = readLahtveeData(1);

costs = cell(length(data.lipidData.metIDs),2);
for i = 1:length(data.lipidData.metIDs)
    %Add exchange reaction for metabolite:
    metID = data.lipidData.metIDs{i};
    model_i = addReaction(model,'r_XXXX','metaboliteList',{metID},'stoichCoeffList',-1);
    
    %Optimize it:
    model_i = changeObjective(model_i,{'r_XXXX'},+1);
    sol     = optimizeCbModel(model_i);
    costs{i,1} = data.lipidData.metAbbrev{i};
    costs{i,2} = num2str(sol.f);
end

rmpath('../data')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
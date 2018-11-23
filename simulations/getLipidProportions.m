%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidProps = getLipidProportions(data,abundances)
%
% Benjamin J. Sanchez. Last update: 2018-11-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lipidProps = getLipidProportions(data,abundances)

molarAbundances = abundances./data.MWs; %(mg/gDW)/(mg/mol) = mol/gDW
lipidProps  = zeros(length(data.allBacks),length(data.allChains));
for i = 1:length(data.metNames)
    pos_b = strcmp(data.allBacks,data.backNames{i});
    for j = 1:length(data.chainNames{i})
        pos_c = strcmp(data.allChains,data.chainNames{i}{j});
        lipidProps(pos_b,pos_c) = lipidProps(pos_b,pos_c) + molarAbundances(i);
    end
end
lipidProps = lipidProps./sum(lipidProps,2)*100;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

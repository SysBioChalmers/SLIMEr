%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L] = sumBioMass(model,comps)
% Calculates breakdown of biomass:
% X -> Biomass fraction without lipids [g/gDW]
% P -> Protein fraction [g/gDW]
% C -> Carbohydrate fraction [g/gDW]
% R -> RNA fraction [g/gDW]
% D -> DNA fraction [g/gDW]
% L -> Lipid fraction [g/gDW]
%
% Benjamin J. Sanchez. Last update: 2018-09-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model,comps)

%Get main fractions:
[P,X] = getFraction(model,comps,'P',0);
[C,X] = getFraction(model,comps,'C',X);
[R,X] = getFraction(model,comps,'R',X);
[D,X] = getFraction(model,comps,'D',X);
[L,X] = getFraction(model,comps,'L',X);

%Add up any remaining components:
bioPos = strcmp(model.rxns,'r_4041');
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X         = X + abundance;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,comps,compType,X)

%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];
rxnName = strrep(rxnName,'P','protein');
rxnName = strrep(rxnName,'C','carbohydrate');
rxnName = strrep(rxnName,'R','RNA');
rxnName = strrep(rxnName,'D','DNA');
rxnName = strrep(rxnName,'L','lipid');

%Add up fraction:
if contains(rxnName,'lipid')
    fractionPos = strcmp(model.rxnNames,[rxnName ' - backbone']);
    subs        = model.S(:,fractionPos) < 0;        %substrates in pseudo-rxn
    F           = -sum(model.S(subs,fractionPos));   %g/gDW
else
    fractionPos = strcmp(model.rxnNames,rxnName);
    comps = comps(strcmp(comps(:,3),compType),:);
    F = 0;
    %Add up all components:
    for i = 1:length(model.mets)
        pos = strcmp(comps(:,1),model.mets{i});
        if sum(pos) == 1
            abundance = -model.S(i,fractionPos)*(comps{pos,2}-18)/1000;
            F         = F + abundance;
        end
    end
end
X = X + F;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

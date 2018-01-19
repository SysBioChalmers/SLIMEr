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
% Benjamín J. Sánchez. Last update: 2018-01-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L] = sumBioMass(model,comps)

bioPos = strcmp(model.rxns,'r_4041');
X      = 0;
P      = 0;
C      = 0;
R      = 0;
D      = 0;
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        abundance = -model.S(i,bioPos)*comps{pos,2}/1000;
        X = X + abundance;
        if strcmp(comps{pos,3},'P')
            P = P + abundance;
        elseif strcmp(comps{pos,3},'C')
            C = C + abundance;
        elseif strcmp(comps{pos,3},'R')
            R = R + abundance;
        elseif strcmp(comps{pos,3},'D')
            D = D + abundance;
        end
    end
end

lipidPos = strcmp(model.rxnNames,'lipid pseudoreaction - backbone');
if sum(lipidPos) == 0
    lipidPos = strcmp(model.rxnNames,'lipid pseudoreaction');
end
subs  = model.S(:,lipidPos) < 0;        %substrates in lipid pseudo-rxn
L     = -sum(model.S(subs,lipidPos));   %lipid composition
X     = X + L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
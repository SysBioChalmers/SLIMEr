%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [model,GAMpol] = changeOtherComp(model,data)
%
% Benjamín J. Sánchez. Last update: 2018-01-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,GAMpol] = changeOtherComp(model,data)

otherData = data.otherData;
fluxData  = data.fluxData;

%Components of biomass:
%        id             MW [g/mol]  class     name
comps = {'s_0404[c]'	89.09       'P'     % A     Alanine         ala
         's_0542[c]'    121.16      'P'     % C     Cysteine        cys
         's_0432[c]'    133.11      'P'     % D     Aspartic acid   asp
         's_0748[c]'    147.13      'P'     % E     Glutamic acid   glu
         's_1314[c]'    165.19      'P'     % F     Phenylalanine   phe
         's_0757[c]'    75.07       'P'     % G     Glycine         gly
         's_0832[c]'    155.15      'P'     % H     Histidine       his
         's_0847[c]'    131.17      'P'     % I     Isoleucine      ile
         's_1099[c]'    146.19      'P'     % K     Lysine          lys
         's_1077[c]'    131.17      'P'     % L     Leucine         leu
         's_1148[c]'    149.21      'P'     % M     Methionine      met
         's_0430[c]'    132.12      'P'     % N     Asparagine      asn
         's_1379[c]'    115.13      'P'     % P     Proline         pro
         's_0747[c]'    146.14      'P'     % Q     Glutamine       gln
         's_0428[c]'    174.2       'P'     % R     Arginine        arg
         's_1428[c]'    105.09      'P'     % S     Serine          ser
         's_1491[c]'    119.12      'P'     % T     Threonine       thr
         's_1561[c]'    117.15      'P'     % V     Valine          val
         's_1527[c]'    204.23      'P'     % W     Tryptophan      trp
         's_1533[c]'    181.19      'P'     % Y     Tyrosine        tyr
         's_0001[ce]'	180.16      'C'     % (1->3)-beta-D-glucan
         's_0004[ce]'	180.16      'C'     % (1->6)-beta-D-glucan
         's_0509[c]'    221.21      'C'     % chitin
         's_0773[c]'    180.16      'C'     % glycogen
         's_1107[c]'    180.16      'C'     % mannan
         's_1520[c]'    342.296 	'C'     % trehalose
         's_0423[c]'    347.22      'R'     % AMP
         's_0526[c]'    323.2       'R'     % CMP
         's_0782[c]'    363.22      'R'     % GMP
         's_1545[c]'    324.18      'R'     % UMP
         's_0584[c]'    331.22      'D'     % dAMP
         's_0589[c]'    307.2       'D'     % dCMP
         's_0615[c]'    345.21      'D'     % dGMP
         's_0649[c]'    322.21      'D'     % dTMP
         's_3714[c]'    852.83      'N'     % heme a
         's_1405[c]'    376.36      'N'     % riboflavin
         's_1467[c]'    96.06       'N'};   % sulphate

%Change given abundances in model:
[~,P,~,R,~] = sumBioMass(model,comps);
bioPos      = strcmp(model.rxns,'r_4041');
for i = 1:length(otherData.metIDs)
    if strcmp(otherData.metIDs{i},'protein')
        fP   = otherData.abundance(i)/P;                            %ratio to scale
        isAA = ~cellfun(@isempty,strfind(model.metNames,'tRNA'));   %protein components
        model.S(isAA,bioPos) = full(model.S(isAA,bioPos))*fP;
        
    elseif strcmp(otherData.metIDs{i},'RNA')
        fR    = otherData.abundance(i)/R;           %ratio to scale
        nucs  = comps(strcmp(comps(:,3),'R'),1);    %RNA components
        for j = 1:length(nucs)
            modelPos = strcmp(model.mets,nucs{j});
            model.S(modelPos,bioPos) = full(model.S(modelPos,bioPos))*fR;
        end
        
    else
        modelPos = strcmp(model.mets,otherData.metIDs{i});
        compPos  = strcmp(comps(:,1),otherData.metIDs{i});
        model.S(modelPos,bioPos) = -otherData.abundance(i)/comps{compPos,2}*1000;
    end
end

%Compute new biomass and lipid fraction:
[X,~,C,~,~] = sumBioMass(model,comps);
lipidPos    = strcmp(model.rxnNames,'lipid pseudoreaction - backbone');
if sum(lipidPos) == 0
    lipidPos = strcmp(model.rxnNames,'lipid pseudoreaction');
end
subs  = model.S(:,lipidPos) < 0;        %substrates in lipid pseudo-rxn
L     = -sum(model.S(subs,lipidPos));   %lipid composition
delta = (X+L)-1;                        %difference to balance

%Balance out mass with all sugars:
mets     = comps(strcmp(comps(:,3),'C'),1);
massPre  = C;
massPost = massPre - delta;
fC       = massPost/massPre;
for i = 1:length(mets)
    modelPos = strcmp(model.mets,mets{i});
    model.S(modelPos,bioPos) = model.S(modelPos,bioPos)*fC;
end

%Estimate maintenance belonging to polymerization:
[~,P,C,R,D] = sumBioMass(model,comps);
GAMpol      = P*37.7 + C*12.8 + R*26.0 + D*26.0;    %Förster 2003 (sup table 8)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
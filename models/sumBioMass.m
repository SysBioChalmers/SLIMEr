%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X = sumBioMass(model)
%
% Benjamín J. Sánchez. Last update: 2018-01-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = sumBioMass(model)

%Components of biomass:
%        id             MW [g/mol]    name
comps = {'s_0404[c]'	89.09       % A     Alanine         ala
         's_0542[c]'    121.16      % C     Cysteine        cys
         's_0432[c]'    133.11      % D     Aspartic acid   asp
         's_0748[c]'    147.13      % E     Glutamic acid   glu
         's_1314[c]'    165.19      % F     Phenylalanine   phe
         's_0757[c]'    75.07       % G     Glycine         gly
         's_0832[c]'    155.15      % H     Histidine       his
         's_0847[c]'    131.17      % I     Isoleucine      ile
         's_1099[c]'    146.19      % K     Lysine          lys
         's_1077[c]'    131.17      % L     Leucine         leu
         's_1148[c]'    149.21      % M     Methionine      met
         's_0430[c]'    132.12      % N     Asparagine      asn
         's_1379[c]'    115.13      % P     Proline         pro
         's_0747[c]'    146.14      % Q     Glutamine       gln
         's_0428[c]'    174.2       % R     Arginine        arg
         's_1428[c]'    105.09      % S     Serine          ser
         's_1491[c]'    119.12      % T     Threonine       thr
         's_1561[c]'    117.15      % V     Valine          val
         's_1527[c]'    204.23      % W     Tryptophan      trp
         's_1533[c]'    181.19      % Y     Tyrosine        tyr
         's_0001[ce]'	180.16      % (1->3)-beta-D-glucan
         's_0004[ce]'	180.16      % (1->6)-beta-D-glucan
         's_0509[c]'    221.21      % chitin
         's_0773[c]'    180.16      % glycogen
         's_1107[c]'    180.16      % mannan
         's_1520[c]'    342.296 	% trehalose
         's_0423[c]'    347.22      % AMP
         's_0526[c]'    323.2       % CMP
         's_0782[c]'    363.22      % GMP
         's_1545[c]'    324.18      % UMP
         's_0584[c]'    331.22      % dAMP
         's_0589[c]'    307.2       % dCMP
         's_0615[c]'    345.21      % dGMP
         's_0649[c]'    322.21      % dTMP
         's_3714[c]'    852.83      % heme a
         's_1405[c]'    376.36      % riboflavin
         's_1467[c]'    96.06};     % sulphate

bioPos = strcmp(model.rxns,'r_4041');
X      = 0;
for i = 1:length(model.mets)
    pos = strcmp(comps(:,1),model.mets{i});
    if sum(pos) == 1
        disp(model.metNames{i})
        disp(num2str(model.S(i,bioPos)*comps{pos,2}/1000))
        X = X + -model.S(i,bioPos)*comps{pos,2}/1000;
    end
end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWs = getMWfromFormula(metFormulas)
%
% Benjamín J. Sánchez. Last update: 2017-11-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MWs = getMWfromFormula(metFormulas)

%Atomic weights:
AWs = {'H' 1.00794
       'C' 12.011
       'N' 14.00674
       'O' 15.9994
       'P' 30.973762
       'S' 32.066};

%Calculate the MW for each formula:
MWs = zeros(size(metFormulas));
for i = 1:length(MWs)
    formula = metFormulas{i};
    elePos  = [regexp(formula,'[A-Z]') length(formula)+1];
    %Find stochiometry of each element and multiply by atomic weight:
    for j = 1:length(elePos)-1
        element = formula(elePos(j));
        AW      = AWs{strcmp(AWs(:,1),element),2}/1000; %[g/mmol]
        stoich  = str2double(formula(elePos(j)+1:elePos(j+1)-1));
        if isnan(stoich)
            stoich = 1;
        end
        MWs(i) = MWs(i) + stoich*AW;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
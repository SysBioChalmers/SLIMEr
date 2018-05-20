%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MWs = getMWfromFormula(metFormulas)
% Returns the MW (in g/mmol) of a chemical formula
%
% Benjamín J. Sánchez. Last update: 2018-05-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MWs = getMWfromFormula(metFormulas)

%Atomic weights:
elements = {'C',   'H',      'O',     'N',      'P',       'S'};
AWs      = [12.011, 1.00794, 15.9994, 14.00674, 30.973762, 32.066];
       
%Calculate the MW for each formula:
MWs = zeros(size(metFormulas));
for i = 1:length(elements)
    %Find stochiometry of each element and multiply by atomic weight:
    stoich = getStoichFromFormula(metFormulas,elements{i});
    MWs    = MWs + stoich*AWs(i)/1000;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
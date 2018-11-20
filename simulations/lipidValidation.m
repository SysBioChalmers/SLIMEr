%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamin J. Sanchez. Last update: 2018-11-18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
close all
addpath('../data')
addpath('../models')

conds  = {'BY4741 wild type - 24ºC','BY4741 wild type - 37ºC', ...
          'BY4741 \itelo1\rm\Delta - 24ºC','BY4741 \itelo1\rm\Delta - 37ºC', ...
          'BY4741 \itelo2\rm\Delta - 24ºC','BY4741 \itelo2\rm\Delta - 37ºC', ...
          'BY4741 \itelo3\rm\Delta - 24ºC','BY4741 \itelo3\rm\Delta - 37ºC'};
Ncond  = length(conds);
errors = zeros(2,Ncond);
Nsim   = 10000;
data   = readEjsingData(1);
data   = convertEjsingData(data,model,false);
lipids = zeros(Ncond*(2*Nsim+1),length(data.metNames)-1);
for i = 1:Ncond
    %Read experimental data for each condition:
    data = readEjsingData(i);
    data = convertEjsingData(data,model,false);
    data.metNames  = data.metNames(1:end-1);        %filter out ergosterol
    data.abundance = data.abundance(1:end-1)*1000;	%mg/gDW
    data.std       = data.std(1:end-1)*1000;        %mg/gDW
    data.MWs       = data.MWs(1:end-1)*1000;        %mg/mol
    
    %Get backbone & chain information:
    data.backNames  = cell(size(data.metNames));
    data.chainNames = cell(size(data.metNames));
    data.allChains  = [];
    for j = 1:length(data.metNames)
        parts = strsplit(data.metNames{j},' ');
        data.backNames{j}  = parts{1};
        data.chainNames{j} = strsplit(parts{2},'-');
        data.chainNames{j} = regexprep(data.chainNames{j},';[0-9]','');
        data.allChains     = [data.allChains data.chainNames{j}];
    end
    data.groups    = findgroups(data.backNames);
    data.allBacks  = unique(data.backNames);
    data.allChains = unique(data.allChains);
    
    %Perform random sampling:
    abundance_modC{i} = lipidRandomSampling(model_corrComp_val{i},data,Nsim);
    abundance_modS{i} = lipidRandomSampling(model_SLIMEr_val{i},data,Nsim);
    
    %Add experimental and simulated lipids:
    Nexp = 2*Nsim*(i-1)+i;
    lipids(Nexp,:) = data.abundance';
    lipids((Nexp+1):(Nexp+Nsim),:)        = abundance_modC{i}';
    lipids((Nexp+Nsim+1):(Nexp+2*Nsim),:) = abundance_modS{i}';
    
    %Averages - stds - errors:
    means       = [mean(abundance_modC{i},2) mean(abundance_modS{i},2)];
    stds        = [std(abundance_modC{i},[],2) std(abundance_modS{i},[],2)];
    diffs       = means - ones(size(means)).*data.abundance;
    errors(:,i) = sum(abs(diffs'),2)/length(diffs);
    
    %Reshape model results:
    x = repmat((1:length(data.metNames))',Nsim,1);
    y = reshape(abundance_modS{i},[length(data.metNames)*Nsim,1]);
    
    %Plot data:
    ymax  = ceil(max(y)/5)*5;
    trans = 0.006;
    if i == 1
        %Get molar proportions:
        propsExp    = getLipidProportions(data,data.abundance);
        propsPerm   = getLipidProportions(data,means(:,1));
        propsSLIMEr = getLipidProportions(data,means(:,2));
        
        %Filter out any unmeasured lipid class:
        names       = data.allBacks(~isnan(sum(propsExp,2)));
        propsPerm   = propsPerm(~isnan(sum(propsExp,2)),:);
        propsSLIMEr = propsSLIMEr(~isnan(sum(propsExp,2)),:);
        propsExp    = propsExp(~isnan(sum(propsExp,2)),:);
        
        %Estimate molar proportions for a hypothetical restrictive model:
        dataCon    = readEjsingData(i);
        dataCon    = convertEjsingData(dataCon,model,true);
        chainMWs   = getMWfromFormula(dataCon.chainData.formulas);	% g/mmol
        chainAbund = dataCon.chainData.abundance;                   % g/gDW
        chainAbund = chainAbund./chainMWs;                          % mmol/gDW
        chainProps = chainAbund/sum(chainAbund)*100;                % mol%
        propsRest  = ones(size(propsExp)).*chainProps';
        
        %Fig S2: Experimental chain distribution
        xlength = 1350;
        figure('position', [100,100,xlength,400])
        color = sampleCVDmap(6);
        barPlot(propsExp,names,'[molar %]',color,100,xlength,[],true);
        legend(data.allChains,'location','eastoutside')
        legend('boxoff')
        
        %Fig 3A: Random sampling at reference conditions
        xlength = 1400;
        figure('position', [100,300,xlength,400])
        barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
        hold on
        scatter(x,y,10,'b','MarkerEdgeAlpha',trans)
        hold off
        
        %Fig S4B: Cumulative distribution of variability
        figure('position', [100,100,500,500])
        stdmax = ceil(max(max(stds)));
        color  = {'y';'b'};
        hold on
        for j = 1:2
            sorted = sort(stds(:,j));
            cumdis = (0:length(sorted))./length(sorted);
            medval = median(stds(:,j));
            h(j) = plot([0;sorted],cumdis,['-' color{j}],'LineWidth',3);
        end
        plotOptions([0 stdmax],[0 1],'Standard Deviation [mg/gDW]', ...
            'Cumulative distribution',[],[],[],[],12)
        modelNames = {'Permissive model';'Enhanced model'};
        legend(h,['\color[rgb]{1 1 0} ' modelNames{1}], ...
                 ['\color[rgb]{0 0 1} ' modelNames{2}],'Location','southeast');
        legend('boxoff')
        hold off
        
        %Fig S4C: Variability breakdown
        figure('position', [100,100,1200,600])
        xticks = (1:length(stds))';
        hold on
        for j = 1:2
            h(j) = plot(xticks,stds(:,j),['o' color{j}],'MarkerFaceColor',color{j});
            %Compute average stds for each backbone group:
            meanstds = splitapply(@mean,stds(:,j),data.groups);
            for k = 1:length(meanstds)
                meanx = xticks(data.groups == k);
                meany = meanstds(k)*ones(size(meanx));
                plot(meanx,meany,['--' color{j}])
                plot([max(meanx)+0.5,max(meanx)+0.5],[0,stdmax],'--k')
            end
            totalmean = mean(stds(:,j));
            plot([0 length(xticks)+1],[totalmean,totalmean],['--' color{j}],'LineWidth',2.5)
        end
        plotOptions([0 length(xticks)+1],[0 stdmax],[],'Standard Deviation [mg/gDW]', ...
            xticks,[],data.metNames,[],6)
        xtickangle(90)
        legend(h,['\color[rgb]{1 1 0} ' modelNames{1}], ...
                 ['\color[rgb]{0 0 1} ' modelNames{2}],'Location','northwest');
        legend('boxoff')
        hold off
        
        %Fig S5: Comparing molar proportions
        figure('position', [100,100,xlength,1000])
        color  = sampleCVDmap(4);
        Nchain = length(propsExp(1,:)); 
        for j = 1:Nchain
            subplot(Nchain,1,j)
            nameChain = ['C' data.allChains{j} ' [%]'];
            props     = [propsRest(:,j) propsPerm(:,j) propsSLIMEr(:,j) propsExp(:,j)];
            barPlot(props,names,nameChain,color,100,xlength,[],false);
            legend('Restrictive Model','Permissive Model','SLIMEr Model', ...
                   'Experimental Data','location','eastoutside')
            legend('boxoff')
        end
        
        figure('position', [100,100,xlength,600])
    else
        xlength = 1000;
    end
    %Fig S6: Random sampling at all conditions
    subplot(Ncond,1,Ncond+1-i)
    barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
    hold on
    scatter(x,y,10,'b','MarkerEdgeAlpha',trans)
    text(80,ymax-5,conds{i})
    hold off
end

%Fig 3B: PCA @ reference conditions
lipids(lipids < 1e-6) = 1e-6;
lipids = log10(lipids);
[loadings,scores,~,~,explained] = pca(lipids);
figure('position', [100,100,600,600])
hold on
rng default
for j = 1:5:Nsim
    condIndexes = randperm(Ncond);
    for i = condIndexes
        intensity = (1/10+i/Ncond*9/10);
        Nexp  = 2*Nsim*(i-1)+i;
        h1(i) = plot(scores(Nexp+(j:j+4),1), scores(Nexp+(j:j+4),2),'o', ...
                     'Color',[1 1 0]*intensity,'MarkerSize',5);
        h2(i) = plot(scores(Nexp+Nsim+(j:j+4),1),scores(Nexp+Nsim+(j:j+4),2),'o', ...
                     'Color',[0 0 1]*intensity,'MarkerSize',5);
    end
    if rem(j,100) == 1
        disp(['Plotting PCA: ' num2str(j) '/' num2str(Nsim) ' simulations'])
    end
end
for i = 1:Ncond
    Nexp  = 2*Nsim*(i-1)+i;
    h3(i) = plot(scores(Nexp,1),scores(Nexp,2),'or', ...
                 'MarkerFaceColor','r','MarkerSize',5,'LineWidth',2);
end
x_lab = ['PC1: ' num2str(explained(1),2) '% of variation'];
y_lab = ['PC2: ' num2str(explained(2),2) '% of variation'];
legend([h1(4),h2(end),h3(end)],'\color[rgb]{0.5 0.5 0} Permissive model', ...
    '\color[rgb]{0   0   1} Enhanced model', ...
    '\color[rgb]{1   0   0} Experimental values','Location','southeast');
legend('boxoff')
plotOptions('','',x_lab,y_lab,'','','','',15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square
hold off

rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

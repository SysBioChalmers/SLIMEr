%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lipidValidation
%
% Benjamín J. Sánchez. Last update: 2018-05-15
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
    data.metNames  = data.metNames(1:end-1);         %filter out ergosterol
    data.abundance = data.abundance(1:end-1)*1000;   %mg/gDW
    data.std       = data.std(1:end-1)*1000;         %mg/gDW
    
    %Perform random sampling:
    abundance_modC{i} = lipidRandomSampling(model_corrComp_val{i},data,Nsim);
    abundance_modS{i} = lipidRandomSampling(model_SLIMEr_val{i},data,Nsim);
    
    %Add experimental and simulated lipids:
    lipids(Nsim*(2*i-2)+i,:) = data.abundance';
    lipids((Nsim*(2*i-2)+i+1):(Nsim*(2*i-1)+i),:) = abundance_modC{i}';
    lipids((Nsim*(2*i-1)+i+1):(Nsim*(2*i)+i),:)   = abundance_modS{i}';
    
    %Averages & errors:
    means       = [mean(abundance_modC{i},2) mean(abundance_modS{i},2)];
    diffs       = means - ones(size(means)).*data.abundance;
    errors(:,i) = sum(abs(diffs'),2)/length(diffs);
    
    %Reshape model results:
    x = repmat((1:length(data.metNames))',Nsim,1);
    y = reshape(abundance_modS{i},[length(data.metNames)*Nsim,1]);
    
    %Plot data:
    ymax  = ceil(max(y)/5)*5;
    trans = 0.006;
    if i == 1        
        %Fig 3: Random sampling at reference conditions
        xlength = 1400;
        figure('position', [100,300,xlength,400])
        barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
        hold on
        scatter(x,y,10,'b','MarkerEdgeAlpha',trans)
        hold off
        figure('position', [100,100,xlength,600])
    else
        xlength = 1000;
    end
    %Fig S3: Random sampling at all conditions
    subplot(Ncond,1,Ncond+1-i)
    barPlot(data.abundance,data.metNames,'[mg/gDW]','r',ymax,xlength,data.std);
    hold on
    scatter(x,y,10,'b','MarkerEdgeAlpha',trans)
    text(80,ymax-5,conds{i})
    hold off
end

%Fig 4: PCA @ reference conditions
lipids(lipids < 1e-6) = 1e-6;
lipids = log10(lipids);
[loadings,scores,~,~,explained] = pca(lipids);
figure('position', [100,100,600,600])
hold on
for i = 1:Ncond
    h1(i) = plot(scores((Nsim*(2*i-2)+i+1):(Nsim*(2*i-1)+i),1), ...
                 scores((Nsim*(2*i-2)+i+1):(Nsim*(2*i-1)+i),2), ...
                 'o','Color',[1 1 0]*(1/10+i/Ncond*9/10),'MarkerSize',5);
    h2(i) = plot(scores((Nsim*(2*i-1)+i+1):(Nsim*(2*i)+i),1), ...
                 scores((Nsim*(2*i-1)+i+1):(Nsim*(2*i)+i),2), ...
                 'o','Color',[0 0 1]*(1/10+i/Ncond*9/10),'MarkerSize',5);
end
for i = 1:Ncond
    h3(i) = plot(scores(Nsim*(2*i-2)+i,1),scores(Nsim*(2*i-2)+i,2), ...
                 'or','MarkerFaceColor','r','MarkerSize',5,'LineWidth',2);
end
size_text = 15;
xlabel(['PC1: ' num2str(explained(1),2) '% of variation'],'FontSize',size_text)
ylabel(['PC2: ' num2str(explained(2),2) '% of variation'],'FontSize',size_text)
legend([h1(4),h2(end),h3(end)],'\color[rgb]{0.5 0.5 0} Permissive model', ...
                               '\color[rgb]{0   0   1} Enhanced model', ...
                               '\color[rgb]{1   0   0} Experimental values', ...
                               'Location','northeast');
legend('boxoff')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'FontSize',size_text)
set(gca,'XColor','k')
set(gca,'YColor','k')
axis square
box on
hold off

rmpath('../data')
rmpath('../models')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
clc;

load('./results/MCMCAll')
MAll = cell(6,1);
MAll{2} = MCMCAll{1};
MAll{3} = MCMCAll{3};
MAll{5} = MCMCAll{2};
MAll{6} = MCMCAll{4};

load('./results/MCMCAll_joint_nodemmock_demHUA')
MAll{1} = MCMCAll{1};
MAll{4} = MCMCAll{2};

Mpar = [6,9,6,10,6,10,6,9,6,10,6,10;... %me1m and me1h for models a-f
    7,10,7,11,7,11,7,10,7,11,7,11;...  %me2m and me2h for models a-f
    8,NaN,8,NaN,8,NaN,8,11,8,12,8,12;...  %me3m and me3h for models a-f
    NaN,11,9,NaN,9,12,NaN,12,9,NaN,9,13];  %dm and dh for models a-f

Mmodel = [1,1,2,2,3,3,4,4,5,5,6,6];

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
for irate = 1:4 %subfigure - rates
    clear MCMCpar
    subplot('Position',[pos(irate),pos(1),xlen,ylen]);
    
    for j = 1:length(Mpar)
        if isnan(Mpar(irate,j)) == 0
            MCMCpar{j} = MAll{Mmodel(j)}.samples(Mpar(irate,j),:)';
            if mod(j,2) == 1
                if isnan(Mpar(irate,j+1)) == 0
                    violin(MCMCpar{j},j,'facecolor',[100,100,100]./255,'edgecolor','none','facealpha',1);
                else
                    violin(MCMCpar{j},j,'facecolor',[255,153,85]./255,'edgecolor','none','facealpha',1);
                end
            else
                violin(MCMCpar{j},j,'facecolor',[135,222,170]./255,'edgecolor','none','facealpha',1);
            end
            hold on
        end
    end

    xlim([0,14])
    ylim([-5,0.5])
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    if irate == 1
        yticks([-5:0])
        ylabel('abundance')
    else
        yticks([])
        set(gca,'ycolor',[1 1 1])
    end
    xticks([]);
end


set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
print('-dpdf',['./figures/ViolinUncertainites'])
DTdata_import;

figure

pos = [0.1,0.32,0.54,0.76];
xlen = 0.4;
ylen = 0.2;

subplot('Position',[pos(1),pos(1),xlen,ylen]);

for imodel = 1
    
    load('./results/MCMCAll_mock_nodem')
    bmedian = median(10.^(MCMCAll{imodel}.samples(1,:)));
    quantileU = quantile(10.^(MCMCAll{imodel}.samples(1,:)),0.75);
    quantileL = quantile(10.^(MCMCAll{imodel}.samples(1,:)),0.25);
    
    for t = 1:40
        Dmedian(t) = 0.5+bmedian*t/(bmedian+t);
        DquantileU(t) = 0.5+quantileU*t/(quantileU+t);
        DquantileL(t) = 0.5+quantileL*t/(quantileL+t);
    end
    
    plot(0+5.5:40+5.5,[0.5,Dmedian],'-','Color','k','Linewidth',1.02)
    hold on
    x = 0+5.5:40+5.5;
    x2 = [x, fliplr(x)];
    inBetween = [[0.5,DquantileL], fliplr([0.5,DquantileU])];
    fill(x2, inBetween,[100,100,100]./255,'LineStyle','none');
    alpha(.3)
    plot(DTdata(1:end-1,2),DTdata(1:end-1,1)/60,'o','Color','k','Markersize',5)
    hold on
end
ylabel('doubling time (hours)')
xlabel('time (h)')
xlim([5.5,45])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
xticks([10,20,30,40])
yticks([0,5,10])

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
print('-dpdf',['./figures/DT_mock_nodem'])

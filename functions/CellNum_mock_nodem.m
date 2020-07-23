%% Figure 2H

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

    t = 0.01:0.01:40;

    Nmedian = Nmedian0.*exp(2.*log(2)./(2*bmedian+1)^2.*(2.*bmedian^2.*log(2.*bmedian.*t+bmedian+t)+2.*bmedian.*t+t-2.*bmedian^2.*log(bmedian)));
    NquantileU = Nmedian0.*exp(2.*log(2)./(2*quantileU+1)^2.*(2.*quantileU^2.*log(2.*quantileU.*t+quantileU+t)+2.*quantileU.*t+t-2.*quantileU^2.*log(quantileU)));
    NquantileL = Nmedian0.*exp(2.*log(2)./(2*quantileL+1)^2.*(2.*quantileL^2.*log(2.*quantileL.*t+quantileL+t)+2.*quantileL.*t+t-2.*quantileL^2.*log(quantileL)));
    
    plot([5.5:0.01:45.5],log10([4096,Nmedian]),'-','Color','k','Linewidth',1.02)
    hold on
    x = [5.5:0.01:45.5];
    x2 = [x, fliplr(x)];
    inBetween = [log10([0.5,NquantileL]), fliplr(log10([0.5,NquantileU]))];
    fill(x2, inBetween,[100,100,100]./255,'LineStyle','none');
    alpha(.3)
%     hold on
%     plot([7,14.75,19.75,27.5,40],log10([40021,107292,115417,221875,292708]),'o','Color','k','Markersize',5)

end
ylabel('log10(number of cells)')
xlabel('time (h)')
xlim([5.5,45])
box off
set(gca,'linewidth',1.02)
set(gca,'FontSize',11)
set(gca,'FontName','Arial')
xticks([10,20,30,40])
ylim([3,6])

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
print('-dpdf',['./figures/CellNum_mock_nodem'])

H4K20_import;
% H4K20dummy_import;

pos = [0.1,0.32,0.54,0.76];
xlen = 0.2;
ylen = 0.2;

figure
for i = 1:4
    subplot('Position',[pos(i),pos(1),xlen,ylen]);
    for j = 1:3
       H4K20DA{i}(:,j) = exp(DA(j).y(:,i)); 
       H4K20DB{i}(:,j) = exp(DB(j).y(:,i));
    end
    MeanH4K20DA{i} = mean(H4K20DA{i}');
    StdH4K20DA{i} = std(H4K20DA{i}');
    MeanH4K20DB{i} = mean(H4K20DB{i}');
    StdH4K20DB{i} = std(H4K20DB{i}');
    
    plot(DA(1).t+5.5,H4K20DA{i},'.','Markersize',15,'Color', [100,100,100]./255,...
        'Linewidth',1.02)
    hold on
    plot(DB(1).t+5.5,H4K20DB{i},'.','Markersize',15,'Color', [135,222,170]./255,...
        'Linewidth',1.02)
    
    xlim([5.5,45])
    ylim([0,1])
%     ylim([0,0.025]) %for inset generation (H4K20me3)
    box off
    set(gca,'linewidth',1.02)
    set(gca,'FontSize',11)
    set(gca,'FontName','Arial')
    xticks([10,20,30,40])
    if i == 1
        yticks([0:0.2:1])
        ylabel('abundance')
    else
        yticks([])
    end
    xlabel('time (h)')
end

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
print('-dpdf',['./figures/FigData'])
% print('-dpdf',['./figures/FigData_Inset'])

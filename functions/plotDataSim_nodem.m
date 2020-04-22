function plotDataSim_nodem

% H4K20_import;
H4K20dummy_import;

load('./parameters/parameters_mock_laplace_mock_MM_1_r1r2r3')
modelsyms1 = 'mock_MM_1_r1r2r3';

model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');

sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
simulateA = eval(sim_name1);

num_cycpar = 1;

r = vertcat([num_cycpar+4+2,num_cycpar+4+2,num_cycpar+4+2]-1,...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+3]-1,...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+2]-1,...
    [num_cycpar+4+2,num_cycpar+4+2,num_cycpar+4+3]-1,...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+4]-1);
modelDefAll = horzcat(repmat(1:num_cycpar+4,5,1),vertcat(r));

M_sub_r = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, 0, 1; 0, 1, 1;...
    1, 1, 1];
M = M_sub_r;

% figure

Models = [4,5];

for imodel = Models
    
    xi = S(imodel).sol.MS.par(:,1);
    modelDef = modelDefAll(imodel,:);
    dist = 'laplace';
    
    n_xi = numel(xi);
    options.llh.distribution = dist;
    options.llh.scale = 'lin';
    options.ami.sensi_meth = 'forward';
    options.ami.atol = 1e-15;
    options.ami.rtol = 1e-8;
    options.ami.sensi = 1;
    options.llh.indA = modelDef;
    
    tsim = linspace(0,40);
    
    D.Y = repmat(DA(1).y(1, :), length(tsim), 1); %data matrix for 'simulate' function (might not be necessary for some users)
    xiA = xi(options.llh.indA);
    solA = simulateA(tsim,xiA,[],D,options.ami);
    
    pos = [0.1,0.32,0.54,0.76];
    xlen = 0.2;
    ylen = 0.2;
    
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
        
        if imodel == 5
            plot(tsim+5.5,exp(solA.y(:,i)),'-','Color','k','Linewidth',1.02)
            hold on
        else
            plot(tsim+5.5,exp(solA.y(:,i)),'-','Color',[100,100,100]./255,'Linewidth',1.02)
            hold on
        end
        
        xlim([5.5,45])
        ylim([0,1])
%         ylim([0,0.025]) %used for inset generation (HK20me3)
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
end

print('-dpdf',['./figures/FigData_mock_MM_1'])
%     print('-dpdf',['./figures/FigData_mock_MM_1_Inset'])


end
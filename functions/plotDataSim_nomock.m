% H4K20_import;
H4K20dummy_import;

% figure

for model = 1:2
    
    load('./parameters/parameters_mockHUA_laplace_mock_MM_1_r1r2r3_HUA_d_r1r2r3') %model e (5)
    
    modelsyms1 = 'mock_MM_1_r1r2r3';
    modelsyms2 = 'HUA_d_r1r2r3';
    
    model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
    amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
    model_syms2 = sprintf('histonesXenopus%s',modelsyms2);
    amiwrap(model_syms2, [model_syms2,'_syms'], './simulation');
    
    sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
    simulateA = eval(sim_name1);
    sim_name2 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,k,D,options)',model_syms2);
    simulateB = eval(sim_name2);
    
    num_cycpar = 1;
    
    M_sub_r = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, 0, 1; 0, 1, 1;...
        1, 1, 1];
    M = horzcat(M_sub_r, ones(8,1));
    
    if model == 1
        Models = 8;
    else
        Models = 5;
    end
    
    for imodel = Models
        
        xi = S(imodel).sol.MS.par(:,1);
        dist = 'laplace';
        
        n_HUA = 6:8;
        count = 0;
        for j = 1:4
            if M(imodel, j) > 0
                count = count + 1;
                n_HUA(j) = 8 + count;
            end
        end
        
        n_xi = numel(xi);
        options.llh.distribution = dist;
        options.llh.scale = 'lin';
        options.ami.sensi_meth = 'forward';
        options.ami.atol = 1e-15;
        options.ami.rtol = 1e-8;
        options.ami.sensi = 1;
        options.llh.indA = 1:8;
        options.llh.indB = [5,n_HUA];
        
        tsim = linspace(0,40);
        
        options.ami.x0 = [];
        options.ami.sx0 = [];
        D.Y = repmat(DA(1).y(1, :), length(tsim), 1); %data matrix for 'simulate' function (might not be necessary for some users)
        Dsim1.Y=DA(1).y(1,:); %artificial data matrix // needed for windows use
        
        xiA = xi(options.llh.indA);
        solA = simulateA(tsim,xiA,[],D,options.ami);
        
        solC = simulateA(5.5,xiA,[],Dsim1,options.ami);
        xiB = xi(options.llh.indB); %parameter assigenment
        options.ami.x0=transpose((solC.x)); %initial values for HUA
        if options.ami.sensi == 1
            options.ami.sx0=[squeeze(solC.sx(:,:,5:8)),zeros(4,1)];
        end
        solB = simulateB(tsim,xiB,solC.x,D,options.ami); %forward simulation
        
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
            
            if imodel == 13
                plot(DA(1).t+5.5,H4K20DA{i},'.','Markersize',15,'Color', [100,100,100]./255,...
                    'Linewidth',1.02)
                hold on
                plot(DB(1).t+5.5,H4K20DB{i},'.','Markersize',15,'Color', [135,222,170]./255,...
                    'Linewidth',1.02)
            end
            hold on
            if imodel == 5
                plot(tsim+5.5,exp(solA.y(:,i)),'-','Color',[100,100,100]./255,'Linewidth',1.02)
                hold on
                plot(tsim+11,exp(solB.y(:,i)),'-','Color',[135,222,170]./255,'Linewidth',1.02)
                hold on
            else
                plot(tsim+5.5,exp(solA.y(:,i)),'-','Color',[200,200,200]./255,'Linewidth',1.02)
                hold on
                plot(tsim+11,exp(solB.y(:,i)),'-','Color',[215,244,227]./255,'Linewidth',1.02)
                hold on
            end
            
            xlim([5.5,45])
            
            ylim([0,1])
%                     ylim([0,0.025]) %used for inset (H4K30me)
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
    
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[2 4 18 15])
    
    print('-dpdf',['./figures/FigData_joint'])
%         print('-dpdf',['./figures/FigData_joint_Inset'])
    
end
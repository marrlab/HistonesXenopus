function plotDataSim(datatype)

%datatype - either 'mock' or 'HUA' for single models or 'mockHUA' for joint
%models 
    
H4K20_import;
% H4K20dummy_import;

switch datatype
    case 'mock'
        load('./parameters/parameters_mock_laplace_mock_MM_1_d1d2d3_r1r2r3')
        
        modelsyms1 = 'mock_MM_1_d1d2d3_r1r2r3';
        
        model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
        amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
        
        sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
        simulateA = eval(sim_name1);
        
        num_cycpar = 1;
        
    case 'HUA'
        load('./parameters/parameters_HUA_laplace_HUA_d1d2d3_r1r2r3')
        
        modelsyms1 = 'HUA_d1d2d3_r1r2r3';
        
        model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
        amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
        
        sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
        simulateA = eval(sim_name1);
        
        num_cycpar = 0;
        
    case 'mockHUA'
        
        load('./parameters/parameters_mockHUA_laplace_mock_MM_1_d_r1r2r3_HUA_d_r1r2r3')
        
        modelsyms1 = 'mock_MM_1_d_r1r2r3';
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
end

d = repmat(num_cycpar+4+1,5,3); %definition for parameters d1_h d2_h d3_h, models 1-5
dd1 = repmat([num_cycpar+4+1,num_cycpar+4+2,num_cycpar+4+2],5,1); %definition for parameters d1_h d2_h d3_h, models 6-10
dd2 = repmat([num_cycpar+4+1,num_cycpar+4+2,num_cycpar+4+1],5,1); %definition for parameters d1_h d2_h d3_h, models 11-15
dd3 = repmat([num_cycpar+4+1,num_cycpar+4+1,num_cycpar+4+2],5,1); %definition for parameters d1_h d2_h d3_h, models 16-20
d1d2d3 = repmat([num_cycpar+4+1,num_cycpar+4+2,num_cycpar+4+3],5,1); %definition for parameters d1_h d2_h d3_h, models 21-25
%definition for parameters r1_h r2_h r3_h, models 1-5
r = vertcat([num_cycpar+4+2,num_cycpar+4+2,num_cycpar+4+2],...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+3],...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+2],...
    [num_cycpar+4+2,num_cycpar+4+2,num_cycpar+4+3],...
    [num_cycpar+4+2,num_cycpar+4+3,num_cycpar+4+4]);
rr1 = r+1; %definition for parameters r1_h r2_h r3_h, models 6-10
rr2 = r+1; %definition for parameters r1_h r2_h r3_h, models 11-15
rr3 = r+1; %definition for parameters r1_h r2_h r3_h, models 16-20
r1r2r3 = r+2; %definition for parameters r1_h r2_h r3_h, models 21-25
modelDefAll = horzcat(repmat(1:num_cycpar+4,25,1),vertcat(horzcat(d,r),...
    horzcat(dd1,rr1),horzcat(dd2,rr2),horzcat(dd3,rr3),...
    horzcat(d1d2d3,r1r2r3)));

M_sub_r = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, 0, 1; 0, 1, 1;...
    1, 1, 1];
M = horzcat(repmat(M_sub_r, 2, 1),[zeros(8,1); ones(8,1)]);

figure

switch datatype
    case 'mock'
        Models = [14,19,20,24,15,4,5,9,25,10];
    case 'HUA'
        Models = [20, 15, 10, 25, 5];
    case 'mockHUA'
        Models = [13, 8, 16, 5];
end

for imodel = Models
    
    xi = S(imodel).sol.MS.par(:,1);
    modelDef = modelDefAll(imodel,:);
    dist = 'laplace';
    
    if isequal(datatype,'mockHUA')
        n_HUA = 6:9;
        count = 0;
        for j = 1:4
            if M(imodel, j) > 0
                count = count + 1;
                n_HUA(j) = 9 + count;
            end
        end
    end
    
    n_xi = numel(xi);
    options.llh.distribution = dist;
    options.llh.scale = 'lin';
    options.ami.sensi_meth = 'forward';
    options.ami.atol = 1e-15;
    options.ami.rtol = 1e-8;
    options.ami.sensi = 1;
    if ~isequal(datatype,'mockHUA')
        options.llh.indA = modelDef;
    elseif isequal(datatype,'mockHUA')
        options.llh.indA = 1:9;
        options.llh.indB = [5,n_HUA];
    end
    
    tsim = linspace(0,40);
    
    if ~isequal(datatype,'mockHUA')
        D.Y = repmat(DA(1).y(1, :), length(tsim), 1); %data matrix for 'simulate' function (might not be necessary for some users)
        xiA = xi(options.llh.indA);
        solA = simulateA(tsim,xiA,[],D,options.ami);
    else
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
            options.ami.sx0=squeeze(solC.sx(:,:,5:9));
        end
        solB = simulateB(tsim,xiB,solC.x,D,options.ami); %forward simulation
    end
    
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
        
        switch datatype
            case 'mock'
                if imodel == 5
                    plot(DA(1).t+5.5,H4K20DA{i},'.','Markersize',15,'Color', [100,100,100]./255,...
                        'Linewidth',1.02)
                end
                plot(tsim+5.5,exp(solA.y(:,i)),'-','Color',[100,100,100]./255,'Linewidth',1.02)
                hold on
            case 'HUA'
                if imodel == 5
                    plot(DB(1).t+5.5,H4K20DB{i},'.','Markersize',15,'Color', [135,222,170]./255,...
                        'Linewidth',1.02)
                end
                hold on
                if imodel == 5
                    plot(tsim+11,exp(solA.y(:,i)),'-','Color',[0,0,0]./255,'Linewidth',1.02)
                    hold on
                else
                    plot(tsim+11,exp(solA.y(:,i)),'-','Color',[135,222,170]./255,'Linewidth',1.02)
                end
                hold on
            case 'mockHUA'
                if imodel == 13
                    plot(DA(1).t+5.5,H4K20DA{i},'.','Markersize',15,'Color', [100,100,100]./255,...
                        'Linewidth',1.02)
                    hold on
                    plot(DB(1).t+5.5,H4K20DB{i},'.','Markersize',15,'Color', [135,222,170]./255,...
                        'Linewidth',1.02)
                end
                hold on
                plot(tsim+5.5,exp(solA.y(:,i)),'-','Color',[200,200,200]./255,'Linewidth',1.02)
                hold on
                plot(tsim+11,exp(solB.y(:,i)),'-','Color',[215,244,227]./255,'Linewidth',1.02)
                hold on
        end
        xlim([5.5,45])
        
        ylim([0,1])
%         ylim([0,0.025]) % used for inset plot generation (H4K20me3)
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
switch datatype 
    case 'mock'
%         print('-dpdf',['./figures/FigData_mock_MM_1_r1r2r3_d1d2d3'])
%     print('-dpdf',['./figures/FigData_mock_MM_1_r1r2r3_d1d2d3_Inset'])
    case 'HUA'
        print('-dpdf',['./figures/FigData_HUA_r1r2r3_d1d2d3'])
%     print('-dpdf',['./figures/FigData_HUA_r1r2r3_d1d2d3_Inset'])
    case 'mockHUA'
%         print('-dpdf',['./figures/FigData_mock_MM_1_d_r1r2r3_HUA_d_r1r2r3'])
%     print('-dpdf',['./figures/FigData_mock_MM_1_d_r1r2r3_HUA_d_r1r2r3_Inset'])
end
end
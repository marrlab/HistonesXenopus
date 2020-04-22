load('./parameters/parameters_mockHUA_laplace_mock_MM_1_d_r1r2r3_HUA_d_r1r2r3')
modelsyms1 = 'mock_MM_1_d_r1r2r3';
modelsyms2 = 'HUA_d_r1r2r3';
dist = 'laplace';
dem = 'yes';

% H4K20_import;
H4K20dummy_import;

model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
model_syms2 = sprintf('histonesXenopus%s',modelsyms2);
amiwrap(model_syms2, [model_syms2,'_syms'], './simulation');

sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
simulateA = eval(sim_name1);
sim_name2 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,k,D,options)',model_syms2);
simulateB = eval(sim_name2);

DB(1).t = DB(1).t-5.5;

M_sub_r = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 0; 1, 0, 1; 0, 1, 1;...
    1, 1, 1];
M = horzcat(repmat(M_sub_r, 2, 1),[zeros(8,1); ones(8,1)]);

for imodel = [5,8,13,16]
    
    clearvars A
    
    xi = S(imodel).sol;
    
    n_HUA = 6:9;
    count = 0;
    for j = 1:4
        if M(imodel, j) > 0
            count = count + 1;
            n_HUA(j) = 9 + count;
        end
    end
    
    % Building a struct covering all sampling options:
    optionsPesto = PestoOptions();
    optionsPesto.MCMC.nIterations = 1e6;
    optionsPesto.MCMC.mode = optionsPesto.mode;
    
    % PT specific options:
    optionsPesto.MCMC.samplingAlgorithm   = 'PT';
    optionsPesto.MCMC.PT.nTemps           = 5;
    
    % Initialize the chains by making use of the preceeding multi-start local
    % optimization, all of them starting from the same point
    optionsPesto.MCMC.theta0 = S(imodel).sol.MS.par(:,1:5);
    for imat = 1:5
        A(:,:,imat) = inv(squeeze(S(imodel).sol.MS.hessian(:,:,imat)));
    end
    optionsPesto.MCMC.sigma0 = A;
    
    parameters = getParameterSamples(S(imodel).sol, @(xi) logLikelihoodXenopusAll(xi,DA,DB,simulateA,simulateB,n_HUA,dist,dem), optionsPesto);
    savesamples = sprintf('parametersSamples%d',imodel);
    save(savesamples,'parameters')
end
% load('parameters_mock_laplace_mock_MM_1_r1r2r3')
load('./parameters/parameters_mock_laplace_mock_MM_1_r1r2r3')
modelsyms1 = 'mock_MM_1_r1r2r3';
dist = 'laplace';
dem = 'no';

H4K20_import
% H4K20_import_server;
% H4K20dummy_import;

model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');

sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
simulateA = eval(sim_name1);

D = DA;

num_cycpar = 1;

fix = num_cycpar+4;
            
model_parameters = [fix+1 0 0 0;...
    fix+2 fix+1 0 0;...
    fix+1 0 fix+2 0;...
    fix+1 0 0 fix+2;...
    0 fix+1 fix+2 fix+3];

%matrix containing the model specific parameters according to the
%definition of the AMICI model
%columns: parameters K20P_0 K20m1_0 K20m2_0 K20m3_0 noise d1_h d2_h d3_h
%r1_h r2_h r3_h
%rows: model (1-25)
%(modelDef)_ij: parameter j occurs in model i as parameter
%|(modelDef)_ij|
%definition for parameters r1_h r2_h r3_h, models 1-5
r = [fix+1,fix+1,fix+1;...
    fix+1,fix+2,fix+2;...
    fix+1,fix+2,fix+1;...
    fix+1,fix+1,fix+2;...
    fix+1,fix+2,fix+3];
modelDef = horzcat(repmat(1:fix,5,1),r);

for imodel = 5
    
    clearvars A
    
    xi = S(imodel).sol;
    
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
        A(:,:,imat) = inv((S(imodel).sol.MS.hessian(:,:,imat)));
    end
    optionsPesto.MCMC.sigma0 = A;
    
    parameters = getParameterSamples(S(imodel).sol, @(xi) logLikelihoodXenopusAll(xi,D,simulateA,modelDef(imodel,:),dist), optionsPesto);
    savesamples = sprintf('parametersSamples%d_mock_nodem',imodel);
    save(savesamples,'parameters')
end
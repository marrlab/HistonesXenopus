%% main master file for optimization, analysis etc

%requires AMICI (https://github.com/ICB-DCM/AMICI) and 
% PESTO (https://github.com/ICB-DCM/PESTO) 

addpath(genpath(pwd))

runparameterestimation = false;
% runparameterestimation = true;

%for function main_histonesXenopusAll:

%datatype = 'mock' / 'HUA'
%dist = 'normal' / 'laplace'
%modelsyms = 'mock_d1d2d3_r1r2r3'
%num_cycpar = 0 / 1 / 2/ 3
%dem = 'yes' / 'no' / 'nomock' / 'noHUA'
%num_models = 1:25 / 1:16 / 1:5
%server = 'yes'  / 'no'
%modelsyms2 = 'HUA_d_r1r2r3'

%modelsyms1 and modelsyms2 contain information about
%model_#ofpossibleparameters 

% eg. HUA_d_r1r2r3 - HUA model with up to three methylation and one
% demethylation rate constants, no initial states are inferred, initial 
% states given as constants k - only used for HUA modelling in joint model

% eg. HUA_d1d2d3_r1r2r3 - HUA model with up to three methylation and three
% demethylation rate constants

% eg. HUA_r1r2r3_nodem - HUA model with up to three methylation rate constants, no
% demethylation rate constants (no demethylation), no initial states are
% inferred, initial states given as constants k- only used for HUA 
% modelling in joint model

% eg. HUA_r1r2r3 - HUA model with up to three methylation rate constants, no
% demethylation rate constants (no demethylation), initial states are
% inferred - only used as seperate HUA model

% eg. mock_d1d2d3_r1r2r3 - mock model with constant cell cycle function and 
% up to three methylation and three demethylation rate constants

% eg. mock_lin_0_d1d2d3_r1r2r3 - mock model with linear cell cycle function
% without determining the offset and up to three methylation and three 
% demethylation rate constants

% eg. mock_lin_d1d2d3_r1r2r3 - mock model with linear cell cycle function
% determining the offset (setting it to 0.5) and up to three methylation 
% and three demethylation rate constants

% eg. mock_MM_0_d1d2d3_r1r2r3 - mock model with Hill function and Hill 
% coefficient 1 as cell cycle function without further prior kowledge 
% and up to three methylation and three demethylation rate constants

% eg. mock_MM_d1d2d3_r1r2r3 - mock model with Hill function and Hill 
% coefficient 1 as cell cycle function with fixed offset of 0.5 and up to 
% three methylation and three demethylation rate constants

% eg. mock_MM_1_d1d2d3_r1r2r3 - mock model with constrained Hill function 
% and Hill coefficient 1 as cell cycle function with fixed offset of 0.5 
% and up to three methylation and three demethylation rate constants

%for further explanations please see the Methods section 

%% HUA

if runparameterestimation
    
    %Gaussian noise model
    main_histonesXenopusAll('HUA','normal','HUA_d1d2d3_r1r2r3',0,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('HUA','laplace','HUA_d1d2d3_r1r2r3',0,'yes',1:25,false)
    
    %Laplacian noise model / no demethylation
    main_histonesXenopusAll('HUA','laplace','HUA_r1r2r3',0,'no',1:5,false)
    
end

%% mock

if runparameterestimation
    
%constant doubling time function
    
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_r1r2r3',1,'no',1:5,false)
    
%linear doubling time function
    
    %linear without offset = 0.5
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_lin_0_d1d2d3_r1r2r3',2,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_lin_0_d1d2d3_r1r2r3',2,'yes',1:25,false)
    
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_lin_0_r1r2r3',2,'no',1:5,false)
    
%linear with offset = 0.5
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_lin_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_lin_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_lin_r1r2r3',1,'no',1:5,false)
    
%Hill function as doubling time
    
    %Hill function with Hill coefficient 1 without offset
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_MM_0_d1d2d3_r1r2r3',3,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_0_d1d2d3_r1r2r3',3,'yes',1:25,false)
    
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_0_r1r2r3',3,'no',1:5,false)
    
%Hill function with Hill coefficient 1 with offset = 0.5
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_MM_d1d2d3_r1r2r3',2,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_d1d2d3_r1r2r3',2,'yes',1:25,false)
    
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_r1r2r3',2,'no',1:5,false)
    
%Hill function with Hill coefficient 1 with offset = 0.5 and constrained
    %Gaussian noise model
    % main_histonesXenopusAll('mock','normal','mock_MM_1_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_1_d1d2d3_r1r2r3',1,'yes',1:25,false)
 
    %without demethylation
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','mock_MM_1_r1r2r3',1,'no',1:5,false)
    
end

%% HUA data but mock models

if runparameterestimation
    
    %constant doubling time
    %Laplacian noise model
    main_histonesXenopusAll('HUA','laplace','mock_d1d2d3_r1r2r3',1,'yes',1:25,false)
    
    %linear doubling time function - without offset = 0.5
    %Laplacian noise model
    main_histonesXenopusAll('HUA','laplace','mock_lin_0_d1d2d3_r1r2r3',2,'yes',1:25,false)
    
    %Hill function as doubling time - without offset = 0.5
    %Laplacian noise model
    main_histonesXenopusAll('HUA','laplace','mock_MM_0_d1d2d3_r1r2r3',3,'yes',1:25,false)
    
end

%% mock data but HUA model

if runparameterestimation
    
    %Laplacian noise model
    main_histonesXenopusAll('mock','laplace','HUA_d1d2d3_r1r2r3',0,'yes',1:25,false)
    
end

%% mock&HUA (joint model) - Hill function with Hill coefficient 1 and offset = 0.5, underlying model ID 5 for both HUA and mock

if runparameterestimation
    
%joint model with demethylation in both mock and HUA
    %Gaussian noise model
    % main_histonesXenopusAll('mockHUA','normal','mock_MM_1_d_r1r2r3',1,'yes',1:16,false,'HUA_d_r1r2r3')
    
    %Laplacian noise model
    main_histonesXenopusAll('mockHUA','laplace','mock_MM_1_d_r1r2r3',1,'yes',1:16,false,'HUA_d_r1r2r3')
    
    
%joint model without demethylation in both mock and HUA
    
    %Laplacian noise model
    main_histonesXenopusAll('mockHUA','laplace','mock_MM_1_r1r2r3',1,'no',1:8,false,'HUA_r1r2r3_nodem')
    
%joint model without demethylation in mock
    
    %Laplacian noise model
    main_histonesXenopusAll('mockHUA','laplace','mock_MM_1_r1r2r3',1,'nomock',1:8,false,'HUA_d_r1r2r3')
    
%joint model without demethylation in HUA
    
    %Laplacian noise model
    main_histonesXenopusAll('mockHUA','laplace','mock_MM_1_d_r1r2r3',1,'noHUA',1:8,false,'HUA_r1r2r3_nodem')
    
end

%% BICs

addpath(genpath(pwd))

files = dir('./parameters/parameters_*.mat');

for ifiles = 1:length(files)
    
    clearvars S DBIC BIC_sort index logPost
    
    load(files(ifiles).name)
    
    BICval = zeros(length(S),1);
    logPost = zeros(length(S),1);
    
    for k = 1:length(S)
        if isempty(S(k).sol) == 0
            BICval(k,1) = S(k).sol.BIC;
            logPost(k,1) = S(k).sol.MS.logPost(1);
        end
    end
    
    [BIC_sort,index]=sort(BICval,'ascend');
    DBIC = BIC_sort-min(BICval);
    
    SBIC{ifiles}.DBIC = DBIC;
    SBIC{ifiles}.BIC_sort = BIC_sort;
    SBIC{ifiles}.index = index;
    SBIC{ifiles}.logPost = logPost;
    SBIC{ifiles}.name = files(ifiles).name;

end

save('./results/SBIC','SBIC')

%% average doubling times

files = dir('./parameters/parameters_*.mat');

% for ifiles = 1:length(files)
for ifiles = 1:length(files)
    
    clearvars -except ifiles ADT ADT_name files
    
    load(files(ifiles).name)
    
    str = extractBetween(files(ifiles).name,"_","_");
    datatype = str{1};
    str = extractBetween(files(ifiles).name,[datatype,'_'],'_');
    distr = str{1};
    str = extractBetween(files(ifiles).name,[distr,'_'],'.mat');
    modelsyms = str{1};
    
    if contains(modelsyms,'mock')
        
        str = extractBetween(modelsyms,'mock_','d');
        if isempty(str)
            str = extractBetween(modelsyms,'mock_','r');
        end
        
        if length(str{1}) > 7
            str = extractBetween(str{1},'M','r');
            str{1} = ['M',str{1}];
        end
        
        DTfun = str{1};
        
        ADT{ifiles} = AverageDoublingTime(DTfun,files(ifiles).name);
        
    end
    
    ADT_name{ifiles} = files(ifiles).name;
    
end

save('./results/ADT','ADT')
save('./results/ADT_name','ADT_name')

%% MCMC sampling for mock model, joint model with demethylation in mock and HUA, joint model with no demethylation in mock but HUA

runsampling = false;
% runsampling = true;

if runsampling
    
    sampling_mock_nodem %model 5 for cell-cycle duration prediction and cell number prediction (Fig 2 F & G)
    
    sampling_joint %joint models with demethylation in mock and HUA 5, 8, 13, 16 (a, b, c, d) for marginal distributions (Fig 4 E)
    
    sampling_joint_nodemmock_demHUA %joint models with demethylation in only HUA 5, 8 (e, f) for marginal distributions (Fig 4 E)
    
end

%% Process MCMC chain / Uncertainties for joint model mock&HUA

runMCMCprocessing = false;
% runMCMCprocessing = true;

if runMCMCprocessing
    
    MCMC_mock_nodem %model 5 for cell-cycle duration prediction and cell number prediction (Fig 2 F & G)
    
    MCMC_joint %joint models with demethylation in mock and HUA 5, 8, 13, 16 (b, c, e, f) for marginal distributions (Fig 4 E)
    
    MCMC_joint_nodemmock_demHUA %joint models with demethylation in only HUA 5, 8 (a, d) for marginal distributions (Fig 4 E)
    
end

%% Median and interquantile ranges for model parameters (joint model b)

%consider all rate constants (r1m, r2m, r3, d, r1h, r2h, rates numbered 6:11)
%comparison between mock- and HUA specific mono and dimethylation rate
%constants and shared demethylation rate constant (comparison of model b, 
%saved in cell 1 of the MCMCAll structure)
%MCMCAll structure = 1st cell - model b
%MCMCAll structure = 2nd cell - model e
%MCMCAll structure = 3rd cell - model c
%MCMCAll structure = 4th cell - model f

P = 6:11;

for ipar = 1:length(P)
    
    load('MCMCAll')
    bmedian(ipar) = median(10.^(MCMCAll{1}.samples(P(ipar),:)));
    quantileU(ipar) = quantile(10.^(MCMCAll{1}.samples(P(ipar),:)),0.75);
    quantileL(ipar) = quantile(10.^(MCMCAll{1}.samples(P(ipar),:)),0.25);
    
end

%% Median and interquantile ranges for model parameters (joint model f)

%consider mock- and HUA-specific demethylation rates (rates numbered 9 and
%13) for model f, saved in the 4th cell of the MCMCAll structure 
%comparison to shared demethylation rate constant in model b 
%MCMCAll structure = 1st cell - model b
%MCMCAll structure = 2nd cell - model e
%MCMCAll structure = 3rd cell - model c
%MCMCAll structure = 4th cell - model f

P = [9,13];

for ipar = 1:length(P)
    
    load('MCMCAll')
    bmedian(ipar) = median(10.^(MCMCAll{4}.samples(P(ipar),:)));
    quantileU(ipar) = quantile(10.^(MCMCAll{4}.samples(P(ipar),:)),0.75);
    quantileL(ipar) = quantile(10.^(MCMCAll{4}.samples(P(ipar),:)),0.25);
    
end

%% Prediction of cell numbers - mock model without demethylation 5 (at 10 and 40 hpf)

%for mock model without demethylation - only sampled MCMC for one mock
%model (hence imodel = 1)

for imodel = 1
    
    load('MCMCAll_mock_nodem')
    bmedian = median(10.^(MCMCAll{imodel}.samples(1,:)));
    quantileU = quantile(10.^(MCMCAll{imodel}.samples(1,:)),0.75);
    quantileL = quantile(10.^(MCMCAll{imodel}.samples(1,:)),0.25);

    t = 0.01:0.01:40;
    cmedian = (0.5+bmedian.*t./(bmedian+t));
    cquU = (0.5+quantileU.*t./(quantileU+t));
    cquL = (0.5+quantileL.*t./(quantileL+t));
    
    Nmedian0 = 4096;
    for t = 0.01:0.01:40
        if t == 0.01
            Nmedian(100*t) = Nmedian0*exp(log(2)/cmedian(100*t)*0.01);
            NquantileU(100*t) = Nmedian0*exp(log(2)/(0.5+cquU*t/(cquU+t))*0.01);
            NquantileL(100*t) = Nmedian0*exp(log(2)/(0.5+cquL*t/(cquL+t))*0.01);
        else
            Nmedian(single(100*t)) = Nmedian(single(100*t-1))*exp(log(2)/cmedian(single(100*t))*0.01);
            NquantileU(single(100*t)) = NquantileU(single(100*t-1))*exp(log(2)/cquU(single(100*t))*0.01);
            NquantileL(single(100*t)) = NquantileL(single(100*t-1))*exp(log(2)/cquL(single(100*t))*0.01);
        end
    end

end

ind10 = find(0.01:0.01:40 == 10-5.5); 
%time point 10 hpf
Nmedian(ind10) 
NquantileU(ind10) 
NquantileL(ind10)

ind40 = find(0.01:0.01:40 == 40-5.5); 
%time point 40 hpf
Nmedian(ind40) 
NquantileU(ind40) 
NquantileL(ind40) 

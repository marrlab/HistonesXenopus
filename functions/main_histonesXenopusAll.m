function S = main_histonesXenopusAll(datatype,dist,modelsyms1,num_cycpar,dem,num_models,server,modelsyms2)

%datatype = 'mock' / 'HUA' / 'mockHUA' - enter the data set you want to be
%modelled (one in case of 'mock' and 'HUA' or two in case of 'mockHUA')

%dist = 'normal' / 'laplace' - determine the underlying noise distribution
%of the experimental data 

%modelsyms1 = 'mock_d1d2d3_r1r2r3' (this is an example, mock model with 
%constant cell cycle function and allowing up to three methylation and 
%demethylation rate constants, see master script for more) - determine the 
%underlying ODE model - if only one data set is modelled enter here the 
%corresponding ODE model - if both data sets are modelled enter here the 
%mock ODE model

%num_cycpar = 0, 1, 2, 3 - determine how many parameters are required to be estimated 
%for the doubling time function - if HUA model: no cycling is assumed and
%hence also no cycling parameter is needed (= 0) - if mock model: the
%cycling parameters vary from 1 to 3 parameters depending on a constant,
%linearly increasing or a gradually increasing and plateauing cell cycle
%function as well as the prior information assumed  

%dem = 'yes' / 'no' / 'nomock' / 'noHUA'- determine whether demethylation is 
%present in the ODE model or not - this has to match the ODE model
%definition(s) of modelsyms1 and if needed modelsyms2 - 'nomock' describes
%the joint model without demethylation in mock but demethylation in HUA -
%'noHUA' describes the joint model without demethylation in HUA but in mock

%num_models = 1:25 - specifiy which model hypothesis should be run (5, 8,
%16, 25 depending on mock/HUA/mockHUA and whether demethylation is present, 
%see Figure 2B and Figure 4B)

%server = false / true - determine whether optimization is run via server (changed
%paths)

%modelsyms2 = 'HUA_d1d2d3_r1r2r3'(this is an example, HUA model allowing up 
%to three methylation and demethylation rate constants, see master script 
%for more) - determine the underlying ODE model of the HUA model in case 
%of joint modeling - has to match demethylation input

%% add paths and install AMICI
if server
    addpath(genpath('/mnt/znas/icb_zstore01/scratch/users/lea.schuh/XenopusIII'));
else
    addpath(genpath('../HistonesXenopus')) 
    addpath(genpath('./tools/AMICI-master'))
    addpath(genpath('./tools/PESTO-master'))
end

clearvars -except datatype dist modelsyms1 num_cycpar dem num_models modelsyms2 server

clc;

installAMICI 

%% import data
if server
    H4K20_import_server;
else
    % H4K20_import;
    H4K20dummy_import;
end
if isequal(datatype,'HUA')
    D = DB;
elseif isequal(datatype,'mock')
    D = DA;
end

%% open .txt document counting the iteration/model TODO: WHY? WHAT IS THAT?
if isequal(datatype,'mockHUA')
    filename = sprintf('count_%s_%s_%s_%s',datatype,dist,modelsyms1,modelsyms2);
else
    filename = sprintf('count_%s_%s_%s',datatype,dist,modelsyms1);
end
fid = fopen([filename,'.txt'],'w');
fprintf(fid,'start');


%% model simulation compilation
if isequal(datatype,'mockHUA')
    model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
    amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
    model_syms2 = sprintf('histonesXenopus%s',modelsyms2);
    amiwrap(model_syms2, [model_syms2,'_syms'], './simulation');
else
    model_syms1 = sprintf('histonesXenopus%s',modelsyms1);
    amiwrap(model_syms1, [model_syms1,'_syms'], './simulation');
end


%% define simulation file
if isequal(datatype,'mockHUA')
    sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
    simulateA = eval(sim_name1);
    sim_name2 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,k,D,options)',model_syms2);
    simulateB = eval(sim_name2);
else
    sim_name1 = sprintf('@(t,xi,k,D,options) simulate_%s(t,xi,[],D,options)',model_syms1);
    simulateA = eval(sim_name1);
end

%% Match timepoints
% the mock and joint models begin at 5.5 hpf - the data points for both mock 
% and HUA are given with respect to this initiation time 
% the HUA model though, starts only 11 hpf, another 5.5 hours after the
% mock and joint model starts - hence the HUA time points need to be shifted 
% by another 5.5. hours to account for the HUA model start
if isequal(datatype,'HUA')
    D(1).t = D(1).t-5.5;
elseif isequal(datatype,'mockHUA')
    DB(1).t = DB(1).t-5.5;
end

%% Determine the correct parameter combinations for all models - shared / specific / joint parameters
% the model definition (syms) files assume the largest possible models e.g.
% three methylation and three demethylation rate constants. Here we assign
% the parameters specifically to the smaller / nested models.

if isequal(datatype,'mockHUA') % joint model
    %M_sub_r = [r1 r2 r3], '1' signals that the according variable is
    %mock/HUA-specific, '0' signals that the according variable is a joint rate
    %constant
    %required only for HUA model
    M_sub_r = [0, 0, 0; %all methylation rate constant m1, m2, m3 are shared between the mock and HUA models in the joint model
               1, 0, 0; %m1 is assigned as a mock / HUA-specific rate constant, m2, m3 are shared between the mock and HUA models in the joint model
               0, 1, 0; %...
               0, 0, 1; 
               1, 1, 0; 
               1, 0, 1; 
               0, 1, 1;
               1, 1, 1]; %all methylation rate constant m1, m2, m3 are mock / HUA-specific in the joint model
           
    %add single HUA/mock-specific or joint demethylation rate constants
    %'1' signals that the according parameter is mock/HUA-specific, 
    %'0' signals that the according parameter is a joint rate
    %constant in the joint model 
    switch dem
        case 'no' % no demethylation in joint model - the metyhlation rate constant assignment defines all models
        M = M_sub_r;
        case 'nomock' %no demethylation in mock - add the HUA specific demethylation rate here
        M = horzcat(M_sub_r, ones(8,1));
        case 'noHUA' %no demethylation in HUA - do not need demethylation assignment for HUA
        M = M_sub_r;
        case 'yes' %demethylation in both mock and HUA - allow for a shared and specific demethylation rate
        M = horzcat(repmat(M_sub_r, 2, 1),[zeros(8,1); ones(8,1)]);
    end
else
    %define shared and specific model parameters for single mock / HUA models
    %in comparison to the mockHUA case and to account for the nested / 
    %smaller models the actual parameter numbers are assigned here instead 
    %of encoded present or not present as 1 and 0:
    %cell_cycle parameters, me1_0, me2_0, me3_0, noise (4+cell_cycle 
    %parameters) are fix parameters inferred by every model
    %this will give us a matrix with dimensions # of models x # of
    %methylation and demethylation parameters (8), where we have 4 possible
    %methylation rate constants m, m1, m2, m3, where m is a shared and m1,
    %m2, m3 a methylation specific rate constant
    %same is true for demethylation 
    switch dem
        case 'yes'
            %determine the demethylation rate constants first 
            %num_cycpar is the number of parameters needed for the cell
            %cycle function
            %for HUA model num_cycpar = 0
            %for mock models num_cycpar between 1 and 3 depending on cell
            %cycle function 
            
            fix = num_cycpar+4;
            
            model_parameters_d = [repmat([fix+1 0 0 0], 5, 1); ... %only one shared demethylation rate constant 
                repmat([fix+2 fix+1 0 0], 5, 1);... %d1 specfific and one shared demethylation rate constant 
                repmat([fix+1 0 fix+2 0], 5, 1); ... %d2 specfific and one shared demethylation rate constant 
                repmat([fix+1 0 0 fix+2], 5, 1);... %d3 specfific and one shared demethylation rate constant 
                repmat([0 fix+1 fix+2 fix+3], 5, 1)]; %three specfific demethylation rate constant 
            
            %determine now the methylation rate constants in the same way
            %as the methylation rate constants are defined after the
            %demethylation rate constants the parameter numbering here depends
            %on the present demethylation parameters
            %get the parameter number up to now
            for id = 1:5
                numd(id) = max(model_parameters_d(5*id-4,:));
            end
            
            %from there add and number the methylation rate constants 
            for ir = 1:5
                model_parameters_dr{ir} = [numd(ir)+1 0 0 0;... %only one shared methylation rate constant 
                    numd(ir)+2 numd(ir)+1 0 0;... %m1 specfific and one shared demethylation rate constant 
                    numd(ir)+1 0 numd(ir)+2 0;...  %m2 specfific and one shared demethylation rate constant
                    numd(ir)+1 0 0 numd(ir)+2;...  %m3 specfific and one shared demethylation rate constant
                    0 numd(ir)+1 numd(ir)+2 numd(ir)+3]; %three specfific methylation rate constant
            end
            
            %set together the methylation rate constants for all models
            model_parameters_r = cell2mat(model_parameters_dr');
            
            %add together the numbered parameters - first the demethylation
            %rate constants and then the methylation rate constants
            model_parameters = horzcat(model_parameters_d, ...
                model_parameters_r);

            %similar to the definition of the parameters we here define the
            %full model according to the ODE definition
            %here the matrix will be # of models x # of ALL model parameters
            %max(13):
            % cell_cycle parameters (0-3), me1_0, me2_0, me3_0, noise,
            %d1, d2, d3, m1, m2, m3
            
            %again first define demethylatin rate constants in model
            d = repmat(fix+1,5,3); %only one shared demethylation rate constant 
            dd1 = repmat([fix+1,fix+2,fix+2],5,1); %d1 specfific and one shared demethylation rate constant 
            dd2 = repmat([fix+1,fix+2,fix+1],5,1); %d2 specfific and one shared demethylation rate constant 
            dd3 = repmat([fix+1,fix+1,fix+2],5,1); %d3 specfific and one shared demethylation rate constant 
            d1d2d3 = repmat([fix+1,fix+2,fix+3],5,1); %three specfific demethylation rate constant
            
            %then define the corresponding methylationa rate constants
            for ir = 1:5
                r{ir} = [numd(ir)+1 numd(ir)+1 numd(ir)+1;... %only one shared methylation rate constant 
                    numd(ir)+1 numd(ir)+2 numd(ir)+2;... %m1 specfific and one shared demethylation rate constant 
                    numd(ir)+1 numd(ir)+2 numd(ir)+1;...  %m2 specfific and one shared demethylation rate constant
                    numd(ir)+1 numd(ir)+1 numd(ir)+2;...  %m3 specfific and one shared demethylation rate constant
                    numd(ir)+1 numd(ir)+2 numd(ir)+3]; %three specfific methylation rate constant
            end
            
            %set together the demethylation and methylation rate constants
            %with the rest of the model (fix parameters eg. me1_0)
            modelDef = horzcat(repmat(1:fix,25,1),vertcat(d,dd1,dd2,dd3,d1d2d3),...
                cell2mat(r'));
        
        case 'no'
            %determine only the methylation rate constants here
            %num_cycpar is the number of parameters needed for the cell
            %cycle function
            %for HUA model num_cycpar = 0
            %for mock models num_cycpar between 1 and 3 depending on cell
            %cycle function 
            fix = num_cycpar+4;
            
            model_parameters = [fix+1 0 0 0;...
                fix+2 fix+1 0 0;...
                fix+1 0 fix+2 0;...
                fix+1 0 0 fix+2;...
                0 fix+1 fix+2 fix+3];
            
            %similar to the definition of the parameters we here define the
            %full model according to the ODE definition
            %here the matrix will be # of models x # of ALL model parameters
            %max(13):
            % cell_cycle parameters (0-3), me1_0, me2_0, me3_0, noise,
            %m1, m2, m3
            r = [fix+1,fix+1,fix+1;...
                fix+1,fix+2,fix+2;...
                fix+1,fix+2,fix+1;...
                fix+1,fix+1,fix+2;...
                fix+1,fix+2,fix+3];
            modelDef = horzcat(repmat(1:fix,5,1),r);
    end
end

%% Parameter definition and multi-start optimization
for i = num_models
    
    %count the iteration/model
    fid = fopen([filename,'.txt'],'w');
    fprintf(fid,'%5d',i);
    
    %definition of parameters occuring in all models and their min/max values
    parameters.name = {'K20m1_0' 'K20m2_0' 'K20m3_0' 'noise'};
    
    params_cyc = {'a','b','d'};
    parameters.name = [params_cyc{1:num_cycpar}, parameters.name];
    parameters.max = [repmat(10,1,num_cycpar), 2, 2, 2, 0];
    parameters.min = [repmat(-10,1,num_cycpar), -4, -4, -4, -2];
    
    if isequal(datatype,'mockHUA')
        
        switch dem
            
            case 'no'
            
                n_HUA = 6:8;

                parameters_fixmock = {'r1_m', 'r2_m', 'r3_m'};
                parameters.name = [parameters.name, parameters_fixmock];
                parameters.max = [parameters.max, 2, 2, 2];
                parameters.min = [parameters.min, -10, -10, -10];

                % vectors containing the potential additional parameters and their boundaries
                % for model HUA
                % (indicated with '_h')
                par_name_diff = {'r1_h', 'r2_h', 'r3_h'};
                par_max_diff = [2, 2, 2];
                par_min_diff = [-10, -10, -10];
            
            case 'nomock'
            
                n_HUA = 6:8;

                parameters_fixmock = {'r1_m', 'r2_m', 'r3_m'};
                parameters.name = [parameters.name, parameters_fixmock];
                parameters.max = [parameters.max, 2, 2, 2];
                parameters.min = [parameters.min, -10, -10, -10];

                par_name_diff = {'r1_h', 'r2_h', 'r3_h', 'd_h'};
                par_max_diff = [2, 2, 2, 2];
                par_min_diff = [-10, -10, -10, -10];
            
            case 'noHUA'
            
                n_HUA = 6:8;

                parameters_fixmock = {'r1_m', 'r2_m', 'r3_m', 'd_m'};
                parameters.name = [parameters.name, parameters_fixmock];
                parameters.max = [parameters.max, 2, 2, 2, -2];
                parameters.min = [parameters.min, -10, -10, -10, -10];

                par_name_diff = {'r1_h', 'r2_h', 'r3_h'};
                par_max_diff = [2, 2, 2];
                par_min_diff = [-10, -10, -10];
            
            case 'yes'
            
                n_HUA = 6:9;

                parameters_fixmock = {'r1_m', 'r2_m', 'r3_m', 'd_m'};
                parameters.name = [parameters.name, parameters_fixmock];
                parameters.max = [parameters.max, 2, 2, 2, 2];
                parameters.min = [parameters.min, -10, -10, -10, -10];

                % vectors containing the potential additional parameters and their boundaries
                % for model HUA
                % (indicated with '_h')
                par_name_diff = {'r1_h', 'r2_h', 'r3_h', 'd_h'};
                par_max_diff = [2, 2, 2, 2];
                par_min_diff = [-10, -10, -10, -10];
            
        end
        
        % creating the correct parameter vector and its boundaries for each run
        count = 0;
        for j = 1:length(par_name_diff)
            if M(i, j) > 0
                parameters.name = horzcat(parameters.name, par_name_diff(j));
                parameters.min = horzcat(parameters.min, par_min_diff(j));
                parameters.max = horzcat(parameters.max, par_max_diff(j));
                count = count + 1;
                % matrix containing the indices of parameters needed for model HUA for the
                % according run
                % listed according to the simplest parameter vector (thus 1:8 or 9 are always included) and the additional
                % parameters (as defined in M)
                switch dem
                    case 'no'
                        n_HUA(j) = 8 + count;
                    
                    case 'nomock'
                    n_HUA(j) = 8 + count;
                    
                    case 'noHUA'
                    n_HUA(j) = 9 + count;
                    
                    case  'yes'
                    n_HUA(j) = 9 + count;
                end
            end
        end
    else
        
        switch dem 
            
            case 'yes'
                %vector containing the model specific parameters and their min/max values
                parameters_name = {'d' 'd1' 'd2' 'd3' 'r' 'r1' 'r2' 'r3'};
                parameters_max = [2, 2, 2, 2, 2, 2, 2, 2];
                parameters_min = [-10, -10, -10, -10, -10, -10, -10, -10];
            case 'no'
                %vector containing the model specific parameters and their min/max values
                parameters_name = {'r' 'r1' 'r2' 'r3'};
                parameters_max = [2, 2, 2, 2];
                parameters_min = [-10, -10, -10, -10];
        end
        
        %parameter vector of model specific parameters for model i
        k=1;
        switch dem
            case 'yes'
                for j = 1:8
                    if model_parameters(i, j) > 0
                        parametersName_part(k) = parameters_name(j);
                        parametersMin_part(k) = parameters_min(j);
                        parametersMax_Part(k) = parameters_max(j);
                        k=k+1;
                    end
                end
            case 'no'
                for j = 1:4
                    if model_parameters(i, j) > 0
                        parametersName_part(k) = parameters_name(j);
                        parametersMin_part(k) = parameters_min(j);
                        parametersMax_Part(k) = parameters_max(j);
                        k=k+1;
                    end
                end
        end
        
        %arrangement of the model specific parameters and their min/max values
        %according to the definition in mock_parameters
        parNum = model_parameters(i,:);
        parNonZero=parNum(parNum~=0)-length(parameters.name);
        parametersName_part = parametersName_part(parNonZero);
        parametersMin_part = parametersMin_part(parNonZero);
        parametersMax_Part = parametersMax_Part(parNonZero);
        
        %full parameter vector with all parameters arranged according to the AMICI
        %model definition for model i
        parameters.name = horzcat(parameters.name,parametersName_part);
        parameters.min =  horzcat(parameters.min,parametersMin_part);
        parameters.max = horzcat(parameters.max,parametersMax_Part);
    end
    
    parameters.number = length(parameters.name);
    
    %Pesto options
    optionsPesto = PestoOptions();
    if server
        optionsPesto.mode = 'silent';
    end
    
    optionsPesto.n_starts = 100;
    
    % local multi start optimization
    if isequal(datatype,'mockHUA')
        S(i).sol = getMultiStarts(parameters,@(xi)logLikelihoodXenopusAll(xi,DA,DB,simulateA,simulateB,n_HUA,dist,dem),optionsPesto);
        S(i).sol.BIC = log(96)*parameters.number-2*S(i).sol.MS.logPost(1);
        % 96 data points
    else
        S(i).sol = getMultiStarts(parameters,@(xi)logLikelihoodXenopusAll(xi,D,simulateA,modelDef(i,:),dist),optionsPesto);
        S(i).sol.BIC = log(48)*parameters.number-2*S(i).sol.MS.logPost(1);
        % 48 data points
    end
end

%% Save results
if isequal(datatype,'mockHUA')
    savename = sprintf('parameters_%s_%s_%s_%s',datatype,dist,modelsyms1,modelsyms2);
else
    savename = sprintf('parameters_%s_%s_%s',datatype,dist,modelsyms1);
end
save(savename, 'S')

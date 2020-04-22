%% Process MCMC chain / Uncertainties for joint model mock&HUA no demethylation for mock but for HUA 

icount = 1;
for imodel = [5,8]
    
    clearvars -except imodel icount CredibleRanges MCMCAll imodel
    
    %get uncertainties
    loadpar = sprintf('./parameters/parametersSamples%d_joint_nodemmock_demHUA',imodel);
    load(loadpar)
    
    %GEWEKE test to dtermine burn-in phase
    frac1 = 0.1;
    frac2 = 0.5;
    [pz,z] = geweke(parameters.S.par,frac1,frac2); %burn-in phase of 10%
    pz %pz should be greater than 0.05
    
    %reduce MCMC chain by burn-in (default 10%)
    [npar,nsamples] = size(parameters.S.par);
    n1 = nsamples*frac1;
    for ipar = 1:npar
        MCMCBurnIn(ipar,:) = parameters.S.par(ipar,n1+1:end);
    end
    
    %reduce MCMC chain further - thinning - only take every 100s sample
    for ipar = 1:npar
        MCMCBurnInThinned(ipar,:) = MCMCBurnIn(ipar,1:100:end);
    end
    
    MCMCAll{icount}.samples = MCMCBurnInThinned;
    
    for jpar = 1:npar
        M(jpar) = mean(10.^(MCMCBurnInThinned(jpar,:)));
        Dev(jpar) = std(10.^(MCMCBurnInThinned(jpar,:)));
    end
    
    MCMCAll{icount}.mean = M;
    MCMCAll{icount}.std = Dev;
    
    icount = icount + 1;
end

save('MCMCAll_joint_nodemmock_demHUA','MCMCAll')
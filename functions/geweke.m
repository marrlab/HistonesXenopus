function [pz,z] = geweke(MCMCchain,frac1,frac2)

%MCMCchain - enter the MCMC chain usded for GEWEKE test and further processing
%frac1 - eg. 0.1 - 10% of burn in phase - fraction of burn in phase 
%frac2 - e.g. 0.5 - comparison to last 50% of MCMC chain 

[npar,nsamples] = size(MCMCchain);

n1 = nsamples*frac1;
n2 = nsamples*frac2;

for i = 1:npar
    S(i,:) = MCMCchain(i,1:n1);
    L(i,:) = MCMCchain(i,end-n2+1:end);
    
    MS(i) = mean(S(i,:));
    SS(i) = std(S(i,:));
    
    ML(i) = mean(L(i,:));
    SL(i) = std(L(i,:));
    
end

z = (MS-ML)./(sqrt(SS+SL));

pz = 2*(1-normcdf(z));

end
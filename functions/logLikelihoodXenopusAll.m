function varargout = logLikelihoodXenopusAll(xi,varargin)
% This functions calculates the logLikelihood function (and its gradient) for
% either one or two data sets (MOCK and MOCK/HUA).
%
% USAGE:
% [logL] = logLikelihoodXenopus()
% [logL,dlogL] = logLikelihoodXenopus()
% [] = logLikelihoodXenopus(xi,DA,simulateA,options) if only one dataset is
% considered
% [] = logLikelihoodXenopus(xi,DA,DB,simulateA,simulateB,options) if both
% datasets are considered
%
%
% Parameters:
%   xi: parameter vector
%   DA,DB: data structs with fields t (time) and y (observables), e.g.:
%       D(1).t = [0,10,20]; D(1).y = [1;1;1]; (n_t x n_y)
%   for one replicate (D(1)) one observables is measured at 3 time points.
%   simulateA,simulateB: simulation functions (provided by AMICI), e.g.:
%       simulateA = @(t,xi,k,D,options) simulateXenopusMOCK(t,xi,[],[],options);
%   options: with fields .llh for the likelihood options and .ami for the
%   AMICI options
%
% Return values:
% logL: logL-likelihood function (be careful with negative log-likelihood and
% positive log-likelihood)
% dlogL: gradient of logL-likelihood function (1 x n_xi)

n_xi = numel(xi);
options.ami.sensi_meth = 'forward';
options.ami.atol = 1e-15;
options.ami.rtol = 1e-8;
options.ami.sensi = 0;
%options.llh.scale = 'lin'; % data on log or linear scale, mainly needed for
%plotting reasons
% options.llh.distribution = 'normal'; % 'laplace'

if nargin > 5
    DA = varargin{1}; % MOCK
    DB = varargin{2}; % HUA
    simulateA = varargin{3}; % MOCK
    simulateB = varargin{4}; % HUA
    n_HUA = varargin{5}; %HUA parameters
    options.llh.distribution = varargin{6};
    dem = varargin{7};
    switch dem
        case 'no'
            options.llh.indA = 1:8;
        
        case 'nomock'  
            options.llh.indA = 1:8;
        
        case 'noHUA'
            options.llh.indA = 1:9;
        
        case 'yes'
            options.llh.indA = 1:9; %indices of parameters for mock data
%         options.llh.indA = 1:6;
    end
    options.llh.indB = [5, n_HUA];  %indices of parameters for HUA data
%     options.llh.indB = [2, n_HUA]; 
    compareflag = 1;
else
    DA = varargin{1};
    simulateA = varargin{2};
    modelDef = varargin{3}; %model parameters
    options.llh.distribution = varargin{4};
    options.llh.indA = modelDef; %indices of parameters for data set A
    compareflag = 0;
end

D.Y = DA(1).y; %artificial data matrix // needed for windows use
Dsim1.Y=DA(1).y(1,:); %artificial data matrix // needed for windows use

try
    nderiv = nargout-1;
    logL = 0;
    if(isfield(options,'sens_ind'))
        if(nderiv>=1)
            dlogL = zeros(length(options.sens_ind),1);
            options.ami.sensi = 1;
        end
    else
        if(nderiv>=1)
            dlogL = zeros(n_xi,1);
            options.ami.sensi = 1;
        end
    end
    
    %% Simulation
    
    xiA = xi(options.llh.indA); %parameter assigenment
    solA = simulateA(DA(1).t,xiA,[],D,options.ami); %forward simulation
    
    if (solA.status<0)
        varargout{1} =NaN;
        if nderiv>=1
            varargout{2} = nan(n_xi,1);
        end
        return; exit
    end
    
    if compareflag %HUA
        %forward simulation of mock model at 5.5h - HUA incubation time
        solC = simulateA(5.5,xiA,[],Dsim1,options.ami);
        xiB = xi(options.llh.indB); %parameter assigenment
        options.ami.x0=transpose((solC.x)); %initial values for HUA
        switch dem 
            case 'no'
            
                if options.ami.sensi == 1
                    options.ami.sx0=[squeeze(solC.sx(:,:,5:8))]; 
                end
            
            case 'nomock'
            
                if options.ami.sensi == 1
                    options.ami.sx0=[squeeze(solC.sx(:,:,5:8)),zeros(4,1)];
                end
            
            case 'noHUA'
            
            if options.ami.sensi == 1
                options.ami.sx0=[squeeze(solC.sx(:,:,5:8))];
            end
            
            case 'yes'
            if options.ami.sensi == 1
                options.ami.sx0=[squeeze(solC.sx(:,:,5:9))]; %initial sensitivities for HUA
                %                 options.ami.sx0=[squeeze(solC.sx(:,:,2:6))];
            end
        end
        solB = simulateB(DB(1).t,xiB,solC.x,D,options.ami); %forward simulation
        if (solB.status<0)
            varargout{1} =NaN;
            if nderiv>=1
                varargout{2} = nan(n_xi,1);
            end
            return;
        end
    end
    
    %% Likelihood evaluation
    switch options.llh.distribution
        case 'normal'
            for idata = 1:length(DA)
                logL = logL - 0.5*sum(sum(~isnan(DA(idata).y)))*log(2*pi)-...
                    nansum(nansum(log(solA.sigmay(~isnan(DA(idata).y))))) -...
                    0.5*nansum(nansum(((solA.y-DA(idata).y)./solA.sigmay).^2));
            end
            if compareflag == 0 %if compareflag == 1 - calculate sensitivities for all parameters together
                if nderiv>=1
                    [N, Edge] = histcounts(options.llh.indA);
                    ind1 = find(N == 1);
                    ind2 = find(N > 1);
                    for iind1 = ind1
                        ind = find(Edge(iind1)+0.5 == options.llh.indA);
                        solA.sy_new(:,:,options.llh.indA(ind))=solA.sy(:,:,ind);
                        solA.ssigmay_new(:,:,options.llh.indA(ind))=solA.ssigmay(:,:,ind);
                    end
                    for iind2 = ind2
                        ind = find(Edge(iind2)+0.5 == options.llh.indA);
                        solsyind = 0;
                        solsigmayind = 0;
                        for iind = ind
                            solsyind = solsyind + solA.sy(:,:,iind);
                            solsigmayind = solsigmayind + solA.ssigmay(:,:,iind);
                        end
                        solA.sy_new(:,:,iind2) = solsyind;
                        solA.ssigmay_new(:,:,iind2) = solsigmayind;
                    end
                    for idata = 1:length(DA)
                        dlogL = dlogL - squeeze(nansum(nansum(bsxfun(@times, ...
                            (solA.y-DA(idata).y)./(solA.sigmay.^2),solA.sy_new),1),2)) - ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solA.sigmay).*(1 - ...
                            (((solA.y-DA(idata).y).^2)./(solA.sigmay.^2)))),solA.ssigmay_new),1),2));
                    end
                end
            end
            
            if compareflag
                for idata = 1:length(DB)
                    logL = logL - 0.5*sum(sum(~isnan(DB(idata).y)))*log(2*pi)-...
                        nansum(nansum(log(solB.sigmay(~isnan(DB(idata).y))))) -...
                        0.5*nansum(nansum(((solB.y-DB(idata).y)./solB.sigmay).^2));
                end
                if nderiv>=1 
                    for idata = 1:length(DA)
                        dlogL(options.llh.indA) = dlogL(options.llh.indA) - squeeze(nansum(nansum(bsxfun(@times, ...
                            (solA.y-DA(idata).y)./(solA.sigmay.^2),solA.sy),1),2)) - ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solA.sigmay).*(1 - ...
                            (((solA.y-DA(idata).y).^2)./(solA.sigmay.^2)))),solA.ssigmay),1),2));
                    end
                    for idata = 1:length(DB)
                        dlogL(options.llh.indB) = dlogL(options.llh.indB) - squeeze(nansum(nansum(bsxfun(@times, ...
                            (solB.y-DB(idata).y)./(solB.sigmay.^2),solB.sy),1),2)) - ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solB.sigmay).*(1 - ...
                            (((solB.y-DB(idata).y).^2)./(solB.sigmay.^2)))),solB.ssigmay),1),2));
                    end
%                     dlogL(ind2) = dlogL(ind2)/2;

                end
            end
            
        case 'laplace'
            for idata = 1:length(DA)
                logL = logL - sum(sum(~isnan(DA(idata).y).*log(2*solA.sigmay))) -...
                    nansum(nansum((abs(solA.y-DA(idata).y)./solA.sigmay)));
            end
            if compareflag == 0 %if compareflag == 1 - calculate sensitivities for all parameters together
                if nderiv>=1
                    [N, Edge] = histcounts(options.llh.indA);
                    ind1 = find(N == 1);
                    ind2 = find(N > 1);
                    for iind1 = ind1
                        ind = find(Edge(iind1)+0.5 == options.llh.indA);
                        solA.sy_new(:,:,options.llh.indA(ind))=solA.sy(:,:,ind);
                        solA.ssigmay_new(:,:,options.llh.indA(ind))=solA.ssigmay(:,:,ind);
                    end
                    for iind2 = ind2
                        ind = find(Edge(iind2)+0.5 == options.llh.indA);
                        solsyind = 0;
                        solsigmayind = 0;
                        for iind = ind
                            solsyind = solsyind + solA.sy(:,:,iind);
                            solsigmayind = solsigmayind + solA.ssigmay(:,:,iind);
                        end
                        solA.sy_new(:,:,iind2) = solsyind;
                        solA.ssigmay_new(:,:,iind2) = solsigmayind;
                    end
                    for idata =  1:length(DA) %sensitivities
                        dlogL = dlogL + squeeze(nansum(nansum(bsxfun(@times, ...
                            sign(DA(idata).y-solA.y)./(solA.sigmay),solA.sy_new),1),2)) + ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solA.sigmay).*((-1) + ...
                            ((abs(solA.y-DA(idata).y))./(solA.sigmay)))),solA.ssigmay_new),1),2));
                    end
                end
            end
            
            if compareflag
                for idata = 1:length(DB)
                    logL = logL - sum(sum(~isnan(DB(idata).y).*log(2*solB.sigmay))) -...
                        nansum(nansum((abs(solB.y-DB(idata).y)./solB.sigmay)));
                end
                if nderiv>=1
                    for idata = 1:length(DA)
                        dlogL(options.llh.indA) = dlogL(options.llh.indA) + squeeze(nansum(nansum(bsxfun(@times, ...
                            sign(DA(idata).y-solA.y)./(solA.sigmay),solA.sy),1),2)) + ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solA.sigmay).*((-1) + ...
                            ((abs(solA.y-DA(idata).y))./(solA.sigmay)))),solA.ssigmay),1),2));
                    end
                    for idata = 1:length(DB)
                        dlogL(options.llh.indB) = dlogL(options.llh.indB) + squeeze(nansum(nansum(bsxfun(@times, ...
                            sign(DB(idata).y-solB.y)./(solB.sigmay),solB.sy),1),2)) + ...
                            squeeze(nansum(nansum(bsxfun(@times,((1./solB.sigmay).*((-1) + ...
                            ((abs(solB.y-DB(idata).y))./(solB.sigmay)))),solB.ssigmay),1),2));
                    end
                end
            end
    end
    
    %% Output assignment
    varargout{1} = logL;
    if nargout >= 2
        varargout{2} = dlogL;
        if sum(isnan(dlogL) > 0)
            varargout{1} = NaN;
            varargout{2} = nan(n_xi,1);
        end
        if sum(isinf(dlogL) > 0)
            varargout{1} = NaN;
            varargout{2} = nan(n_xi,1);
        end
    end
    
catch
    varargout{1} = NaN;
    if nderiv>=1
        varargout{2} = nan(n_xi,1);
    end
end



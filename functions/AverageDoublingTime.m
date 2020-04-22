function ADT = AverageDoublingTime(DTfun,file)

%DTfun = doubling time function according to which the doubling time is
%calculated ('', 'lin_', 'lin_0_', 'MM_', 'MM_0_', 'MM_1_') - either
%constant, linearly increasing or gradually increasing and plateauing with
%different inserted prior knowledge
%file = parameter file name for loading

load(file) %load the corresponding parameter file with inferred parameters

% TODO: CAN YOU DEFINE UP HERE SOMETHING LIKE
% p1 = S(imodel).sol.MS.par(1,1);
% p2 = ...
% something like t = 0:45-5.5
% AND REPLACE IT THROUGHT THIS FUNCTION? TODOEND

%% constant doubling time function
if isempty(DTfun)
    for imodel = 1:length(S)
        ADT(imodel) = 10.^(S(imodel).sol.MS.par(1,1));
    end
    %% linear doubling time function
elseif isequal(DTfun,'lin_') ||  isequal(DTfun,'lin_0_')
    for imodel = 1:length(S)
        switch DTfun
            case 'lin_0_' % linear doubling time function without offset
                ADT(imodel) = mean(10.^(S(imodel).sol.MS.par(1,1))+10.^(S(imodel).sol.MS.par(2,1))*[45-5.5]);
            case 'lin_' % linear doubling time function with offset = 0.5
                ADT(imodel) = mean(0.5+10.^(S(imodel).sol.MS.par(1,1))*[0:45-5.5]);
        end
    end
    
    %% Hill function with coefficient 1
elseif isequal(DTfun,'MM_') ||  isequal(DTfun,'MM_0_') || isequal(DTfun,'MM_1_')
    str = extractBetween(file,DTfun,'.mat');
    for imodel = 1:length(S)
        switch DTfun
            case 'MM_0_' % without offset
                ADT(imodel) = mean(10.^(S(imodel).sol.MS.par(1,1))+...
                    (10.^(S(imodel).sol.MS.par(2,1))*[0:45-5.5])./...
                    (10.^(S(imodel).sol.MS.par(3,1))+[0:45-5.5]));
            case 'MM_' % offset = 0.5
                ADT(imodel) = mean(0.5+(10.^(S(imodel).sol.MS.par(1,1))*[0:45-5.5])./...
                    (10.^(S(imodel).sol.MS.par(2,1))+[0:45-5.5]));
            case 'MM_1_' % offset = 0.5 and constrained
                if isempty(S(imodel).sol) == 0
                    ADT(imodel) = mean(0.5+(10.^(S(imodel).sol.MS.par(1,1))*[0:45-5.5])./...
                        (10.^(S(imodel).sol.MS.par(1,1))+[0:45-5.5]));
                else
                    ADT(imodel) = 0;
                end
        end
    end
    
end

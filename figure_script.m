%file creating all figures
addpath(genpath(pwd))
%% Figure 1B 

plotData

%% Figure 2E

plotDataSim('mock')
plotDataSim_nodem

%% Figure 2F 

DT_mock_nodem

%% Figure 2G

CellNum_mock_nodem

%% Figure 3C

plotDataSim('HUA')

%% Figure 4E

plotDataSim('mockHUA')
plotDataSim_nomock

%% Figure 4F

plotUncertainty

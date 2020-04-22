# Functions

All functions used in the analysis and for plotting the figures/results can be found here. <br>
A brief explanation of all functions can be found here. For more detailed information pleaase have a look at the documentation in the scripts themselves. <br>
<br>

[acf.m](acf.m) - computes autocorrelations through p lags <br>
[AverageDoublingTime.m](AverageDoublingTime.m) - computes the average cell-cylce durations <br>
[CellNum_mock_nodem.m](CellNum_mock_nodem.m) - computes the predicted cell numbers over time for the mock model with three methylation rate constant and no demethylation (Figure 2G) <br>
[DT_mock_nodem.m](DT_mock_nodem.m) - computes the cell-cycle duration over time for the mock model with three methylation rate constant and no demethylation (Figure 2F) <br>
[geweke.m](geweke.m) - GEWEKE test for input MCMC sampling <br>
[logLikelihoodXenopusAll.m](logLikelihoodXenopusAll.m) - computes the loglikelihood function and its gradients for optimization  <br>
[main_histonesXenopusAll.m](main_histonesXenopusAll.m) - main script for parameter estimation/optimization  <br>
[MCMC_joint_nodemmock_demHUA.m](MCMC_joint_nodemmock_demHUA.m) - coorects MCMC sampling for joint models a and d (Figure 4C) - best performing joint models with demethylation only in HUA (used for creating the violin plots in Figure 4F)  <br>
[MCMC_joint.m](MCMC_joint.m) - corrects MCMC sampling for joint models b, c, e and f (Figure 4C) - best performing joint models with demethylation in both mock and HUA (used for creating the violin plots in Figure 4F)  <br>
[MCMC_mock_node.m](MCMC_mock_nodem.m) - corrects MCMC sampling for mock model with three methylation rate constant and no demethylation - used for computing the average cell-cycle duration in Figure 2F and the prediction of cell numbers in Figure 2G  <br>
[plotData.m](plotData.m) - plots the mock and HUA data (Figure 1B)  <br>
[plotDataSim_nodem.m](plotDataSim_nodem.m) - plots mock data with model fits for mock models with demethylation (Figure 2E)  <br>
[plotDataSim_nomock.m](plotDataSim_nomock.m) - plots mock and HUA data and joint model fits for joint models a and d both with demethlyation only present in HUA and not in mock (Figure 4E)  <br>
[plotDataSim.m](plotDataSim.m) - plots mock and HUA data and model fits for mock, HUA and joint models with demethylation (Figure 2E, Figure 3C and Figure 4E)  <br>
[plotUncertainty.m](plotUncertainty.m) - plots the marginal distributions (MCMC sampling) of the rate constants of the best performing joint models as violin plots in Figure 4F  <br>
[sampling_joint_nodemmock_demHUA.m](sampling_joint_nodemmock_demHUA.m) - computes the MCMC samples of joint models a and d both allowing for demethylation in HUA only and not in mock (used for creating the violin plots in Figure 4F)  <br>
[sampling_joint.m](sampling_joint.m) - computes the MCMC samples of joint models b, c, e and f allowing for demethylation in noth mock and HUA (used for creating the violin plots in Figure 4F)  <br>
[sampling_mock_nodem.m](sampling_mock_nodem.m) - computes the MCMC samples of mock model with three methylation rate constants but no demethylation - used for computing the average cell-cycle duration in Figure 2F and the prediction of cell numbers in Figure 2G  <br>
[violin.m](violin.m) - creates violin plots in in Figure 4F  <br>

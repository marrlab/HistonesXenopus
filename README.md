# HistonesXenopus

contains code and data accompanying 

__Computational modeling reveals cell-cycle dependent kinetics of H4K20 methylation states during Xenopus embryogenesis__

Lea Schuh<sup>1,2</sup>, Carolin Loos<sup>1,2,3,4</sup>, Daniil Pokrovsky<sup>5</sup>, Axel Imhof<sup>5</sup>, Ralph A.W. Rupp<sup>5</sup>, Carsten Marr<sup>1,*</sup>

<sub><sup>
<sup>1</sup>Helmholtz Zentrum München-German Research Center for Environmental Health, Institute of Computational Biology, Neuherberg, 85764, Germany <br>
<sup>2</sup>Department of Mathematics, Technical University of Munich, Garching, 85748, Germany <br>
<sup>3</sup>Ragon Institute of MGH, Harvard, MIT, Cambridge, MA 02139, USA <br>
<sup>4</sup>Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA <br>
<sup>5</sup>Department of Molecular Biology, Ludwig-Maximilians-Universität München, Planegg-Martinsried, 82152, Germany <br>
*carsten.marr@helmholtz-muenchen.de <br>
</sup></sub>

## Summary

DNA replication during cell division leads to dilution of histone modifications and can thus affect chromatin-mediated gene regulation, raising the question of how the cell-cycle shapes the histone modification landscape, particularly during embryogenesis. We tackled this problem by manipulating the cell-cycle during early Xenopus laevis embryogenesis and analyzing in vivo histone H4K20 methylation kinetics. The global distribution of un-, mono- di- and tri-methylated histone H4K20 was measured by mass spectrometry in normal and cell-cycle arrested embryos over time. Using multi-start maximum likelihood optimization and quantitative model selection, we found that three specific biological methylation rate constants were required to explain the measured H4K20 methylation state kinetics. While demethylation is essential for regulating H4K20 methylation kinetics in non-cycling cells, demethylation is very likely dispensable in rapidly dividing cells of early embryos, suggesting that cell-cycle mediated dilution of H4K20 methylation is an essential regulatory component for shaping its epigenetic landscape during early development. <br>

Required software and toolboxes:

- MATLAB (R2017a)
- [AMICI](https://github.com/ICB-DCM/AMICI) version 0.10.4, for model definition and model simulations 
- [PESTO](https://github.com/ICB-DCM/PESTO version 1.1.9 for parameter estimation and MCMC sampling

Both toolboxes, AMICI and PESTO, can be found in the folder [tools](tools) but need to be unzipped. Additionally, AMICI creates '.mex' files and requires a C/C++ compiler. The analysis was performed on a macOS Catalina version 10.15.3. <br>

## Folder Structure

Here you can find the list of folders

- [data](data)
- [figures](figures)
- [images](images)
- [parameters](parameters)
- [results](results)
- [syms](syms)
- [tools](tools)

## How to run the code

There are mainly two important scripts (i) for the analysis summarized in [master_script.m](master_script.m) (see Workflow) and (ii) for the figures summarized in [figure_script.m](figure_script.m) (see Figures).

## Data

The H4K20 methylation proportions of both mock and HUA over time can be found in the file [data/H4K20states.xlsx](data/H4K20states.xlsx). Until the original data set will be published we will provide dummy data in [data/H4K20states_dummy.xlsx](data/H4K20states_dummy.xlsx). Upon publication of the data we will upload the the true data here. The doubling times of cells during Xenopus embryogenesis found in literature are summarized in [data/DTdata.xlsx](data/DTdata.xlsx). The corresponding import scripts can be found in the same folder. The images of mock and HUA embryos corresponding to the different developmental stages (Figure 1A) can be found in [images](images).

## Workflow

The workflow to the whole manuscript can be found in the script [master_script.m](master_script.m).
The model files for the compilation with AMICI can be found in the folder [syms](syms), the compiled model files in the folder [simulation](simulation).

## Figures
[figure_script.m](figure_script.m) recreates all figure panels of the paper. The created figure files are saved to and can be found in the folder [figures](figures).

## Results

All estimated parameter sets of the optimization and the MCMC sampling results of selected models can be found in the folder [parameters](parameters). Summaries of the average doubling times of all mock models and the BIC values of all mock and HUA models as well as the processed MCMC samples can be found in the folder [results](results).

## Supplementary Information

The summary of all model optimizations (rank, delta BIC, BIC, model ID adn potentially the average doubling time) can be found in [HistonesXenopus_results.xlsx](HistonesXenopus_results.xlsx).

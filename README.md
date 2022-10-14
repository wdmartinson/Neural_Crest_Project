# Codes for neural crest cell migration / ECM remodeling model

**Release date:** 14 October 2022

## Overview:
This repository contains codes used to simulate the neural crest agent-based mathematical described in the article "Dynamic Fibronectin Assembly and Remodeling by Leader Neural Crest Cells Prevents Jamming in Collective Cell Migration" by WD Martinson, R McLennan, JM Teddy, MC McKinney, LA Davidson, RE Baker, HM Byrne, PM Kulesa, and PK Maini. The model is implemented in C++, using the PhysiCell library (Ghaffarizadeh et al., 2018; version 1.7.1). Codes for creating violin plots have also been used, and have been downloaded from an existing repository created by Bechtold (2016). Additional codes for statistical analysis with partial rank correlation coefficients were created by the authors but inspired by older implementations by Marino et al. (2008, http://malthus.micro.med.umich.edu/lab/usadata/). Implementations for extended Fourier Amplitude Sensitivity Testing (eFAST) relied partially on the Python library SALib (Herman and Usher, 2017). Some summary statistics rely on the Python topological data analysis library Ripser (https://ripser.scikit-tda.org/en/latest/) to generate persitence diagrams. However, results using these tools were not included in the final paper.

The key codes for generating summary statistics are the main_Latin_Hypercube_Sampling_PRCC_Analysis.m, main_eFAST.m, and the main_puncta_experiments.m file. These generate all of the data files and summary statistics that were used to inform the paper. A Dryad repository (DOI: 10.5061/dryad.69p8cz958) stores the particular data files presented in the article.

The main_consistency_analysis.m file runs a script to determine how many realizations of the agent-based model are needed for each parameter regime in order to adequately characterize the distribution of summary statistics. We followed a procedure outlined in Hamis et al. (2021) to perform this analysis, known as consistency analysis.

**Reference:**
WD Martinson, R McLennan, JM Teddy, MC McKinney, LA Davidson, RE Baker, HM Byrne, PM Kulesa, PK Maini, Dynamic Fibronectin Assembly and Remodeling by Leader Neural Crest Cells Prevents Jamming in Collective Cell Migration, 2022. arXiv:2209.07794 (https://arxiv.org/abs/2209.07794)

A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. DOI: [10.1371/journal.pcbi.1005991](https://dx.doi.org/10.1371/journal.pcbi.1005991)

J Herman and W Usher, SALib: an open-source Python library for sensitivity analysis, J. Open Source Softw. 2(9), 97, 2017. DOI: [10.21105/joss.00097]

S Marino, IB Hogue, CJ Ray, DE Kirschner, A methodology for performing global uncertainty and sensitivity analysis in systems biology, J. Theor. Biol. 254(1): 178-196, 2008. DOI: [10.1016/j.jtbi.2008.04.011]

Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project. https://github.com/bastibe/Violinplot-Matlab, DOI: [10.5281/zenodo.4559847]

S Hamis, S Stratiev, G Powathil, Uncertainty and sensitivity analyses methods for agent-based models: an introductory review, The Physics of Cancer, 1-37, 2021. DOI: [10.1142/9789811223495_0001]

### To download Python dependencies:
Open a terminal on your computer. Type the following commands:
pip3 install numpy scipy pandas matplotlib sklearn xml
pip3 install ripser
pip3 install SALib

### Important information for using the PhysiCell library
Visit http://MathCancer.org/blog and http://physicell.org/tutorials/ for installation instructions, latest tutorials, and help.

OSX users must set up gcc and OpenMP on their system in order to run PhysiCell. In addition, they must define the PHYSICELL_CPP system variable. See http://www.mathcancer.org/blog/setting-up-gcc-openmp-on-osx-homebrew-edition/ for more details.

Windows users must follow similar procedures. Information can be found at http://www.mathcancer.org/blog/setting-up-a-64-bit-gcc-environment-on-windows/.

To compile and run the executable file for simulating one model realization, type the following after installing PhysiCell:
make clean
make reset
make ncc-puncta-experiments
make
./project-NCC-puncta-experiment

Parameter values can be changed by updating the PhysiCell_settings.xml file in the config folder of the parent directory.

Functions for the agent-based model can be changed by changing files in the NCC_puncta_experiments folder. Note that you can update any custom modules as well as the default parameter set that is set upon compilation of the executable here.

### Key makefile rules:

**make**               : compiles the current project. If no
                     project has been defined, it first
                     populates the cancer heterogeneity 2D
                     sample project and compiles it

**make \[project-name\]**: populates the indicated sample project.
                     Use "make" to compile it.

  \[project-name\] choices:
    puncta-model-experiments

**make list-projects** : list all available sample projects

**make clean**         : removes all .o files and the executable, so that the next "make" recompiles the entire project

**make data-cleanup**  : clears out all simulation data

**make reset**         : de-populates the sample project and returns to the original PhysiCell state. Use this when switching to a new PhysiCell sample project.


**Homepage for PhysiCell:**     http://PhysiCell.MathCancer.org

**Downloads for PhysiCell:**    http://PhysiCell.sf.net

**Support for PhysiCell:**      https://sourceforge.net/p/physicell/tickets/

**Tutorials for PhysiCell:**    http://www.mathcancer.org/blog/physicell-tutorials/

**Latest info for PhysiCell:**  follow [@MathCancer](https://twitter.com/MathCancer) on Twitter (http://twitter.com/MathCancer)

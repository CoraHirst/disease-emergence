This `disease-emergence/code` folder contains 4 subfolders: 

  1) "data analysis",
  2) "exploratory",
  3) "functions", and 
  4) "project".


`disease-emergence/code/exploratory` contains information about the model underlying the work presented in Antia et al, Nature, 2003. 

`disease-emergence/code/data analysis` contains a subfolder `congo`, in which there are two .qmd documents.
These documents contain the code required to calculate the waning of the fraction immune to smallpox from demographic data on the DRC from 1980-2023 and project continued waning to 2080. See `disease-emergence/code/data analysis/README.md` for more information. 

`disease-emergence/code/project` contains a subfolder `paper figures` in which there is a single .qmd document.
This document contains the code required to generate all manuscript and supplemental figures. 

`disease-emergence/code/functions` contains three .R files. 
`matrix_model.R` contains a function for the leslie matrix required for modeling the change in population size by age over time in the DRC. 
`analysis_functions.R` contains functions for reading the demographic data of the DRC and for calculating the fraction immune over time.
`functions.R` contains all functions related to calculating the nonextinction probabilities of the multi-type branching process functions described in SI 1. 
Relevant functions are loaded from these sources in the data analysis and generation of paper figures files.

Project workflow: The two .qmd documents in `disease-emergence/code/data analysis` should be run before `disease-emergence/code/project/paper figures/paper_figures.qmd`

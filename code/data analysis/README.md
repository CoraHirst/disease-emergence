This `disease-emergence/code/data analysis` folder contains 1 subfolder: `disease-emergence/code/data analysis/congo`.

There are two .qmd documents contained in this subfolder: 
  1) `waning_population_immunity_1.qmd`,
  2) `matrix_model_projections_2.qmd`, and 

NOTE: These files must be run BEFORE generating paper figures. 


Description: 

`waning_population_immunity_1.qmd` describes the source of age-structured demographic data from the DRC, cleans and formats this data to parameterize the leslie matrix model developed in `matrix_model_projections_2.qmd`. 
functions loaded from source `disease-emergence/code/functions/analysis_functions.R`
Input: Demographic data on the DRC from the WHO is saved to `disease-emergence/data/congo/congo_demographics/1980-2023.xlsx`
Output: Birth and death rates, as well as the initial population conditions in 1980, are saved to `disease-emergence/data/congo/formatted data/`: 
  formatted birth rates by age: `birth_rates.csv`
  formatted death rates by age: `death_rates.csv`
  initial population size by age: `pop_1980.csv`

`matrix_model_projections_2.qmd` describes the leslie matrix model saved in `disease-emergence/code/functions/matrix_model.R`, and runs this model to project the change in number of individuals by age over time in the DRC. From this information, the proportion of individuals in the population with smallpox immunity is calculated.
functions loaded from source `disease-emergence/code/functions/matrix_model.R`
Input: Birth and death rates, as well as the initial population conditions in 1980, are saved to `disease-emergence/data/congo/formatted data/`
Output: Age-stratified population matrix output from 1980-2080, and matrix containing the proprotion of individuals suriving with immunity from 1980-2080, are saved to `disease-emergence/results/data analysis/congo/`:
  population matrix: `disease-emergence/results/data analysis/congo/poulation_matrix.csv` 
  proportion immune through time: `disease-emergence/results/data analysis/congo/proportion_immune.csv` 

The files should be run in the following order: 
`waning_population_immunity_1.qmd` then
`matrix_model_projections_2.qmd`

# AR2-Capture-Heterogeneity-Simulation
Simulation code used for assessing the importance of estimating capture heterogeneity while modelling AR(2) processes. 
Results published in Nicolau, Sørbye & Yoccoz (2020). Ecol Evol DOI: 10.22541/au.158385561.11117337 

### The Simulation Data Analysis folder contains 3 functions scripts:
# "0_1_functions_data_generation.R"
# Functions used to generate capture process counts and capture histories for the different time points

# "0_2_functions_for_data_analysis.R"
# Functions used to analyse the generated time series using the different methods

# "0_3_functions_for_run_sims_inparallel.R"
# Functions that put together the previous scripts into running a full simulation setup and estimate the AR coefficients underlying the abundance time series

### And 1 script which run the full simulation code:
# "1_run_simulations.R"
# Script which runs the full simulation exercise, ran in parallel, as well as a function to summarize the different simulations

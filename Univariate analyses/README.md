# Code

- *Splittail assessment functions.R*: Helper functions for Splittail case study
- *Splittail data processing.R*: Process data for model-fitting, and create the sampling effort reduction scenarios
- *Splittail variance analysis.R*: Variance analysis on Splittail data
- *Splittail model fitting.R*: Fit the full and all reduced models
- *Splittail model processing.R*: Evaluate and process the full and reduced models. Create plots for the results of the Splittail case study
- *Additional plots for technical chapter.R*: Create simulation plots for conceptual diagram and Splittail distribution plot for the technical chapter
- *Offset test.R*: Test whether the models would be improved by incorporating sampling effort via an offset of the log of `Tow_area`. 

# Data

All files can be loaded with the `load` function. 

- *Split data.Rds*: The processed Splittail data used to fit all models
- *Full model local trends.Rds*: Processed full model output with 95% credible intervals
- *Reduced model proportions.Rds*: Processed reduced models, specifying the proportion of posterior draws falling within the 95% credible intervals of the full model

# Figures

- *Distribution plots*: A folder of plots of model predictions from the full and each reduced model.
- *Full_sim.png*: Simulated full model for the conceptual diagram
- *Reduced_sim.png*: Simulated reduced models for the conceptual diagram
- *Reduced_sim_overlay.png*: Simulated reduced model over the simulated full model for the conceptual diagram
- *Splittail 0.5 station cut 1 of 2 ribbon example.png*: Plot of example model posteriors from the full and a reduced model
- *Splittail map.png*: Map of Splittail abundance across all sampling stations
- *Splittail reduced model replicates.png*: Data reduction simulation results for each replicate simulation
- *Splittail reduced model summarized.png*: Data reduction simulation results summarized for each season
- *Variance plot.png*: Parameter estimates and 95% credible intervals from the variance analysis
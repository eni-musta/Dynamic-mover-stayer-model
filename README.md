# Dynamic discrete-time mover-stayer-model

Mover-stayer models provide a statistical framework for analysing heterogeneous populations in which some individuals exhibit persistent behaviour ("stayers") while others transition between states over time ("movers"). They are particularly used to understand population dynamics in fields such as social sciences and economics. However, existing mover-stayer models assume a fixed stayer fraction over time (i.e., the mover-stayer status is determined at baseline) and can only incorporate time-dependent covariates in the transition probabilities for the movers.  In this paper we present a novel dynamic version of the mover-stayer model for discrete time-to-event data, allowing individuals who are initially susceptible (potential movers) to become stayers over time based on evolving circumstances. Our approach incorporates both time-fixed and exogenous time-varying covariates to model transition probabilities between the states of potential mover, mover and stayer, through a multinomial logistic model. Estimation is performed using the maximum likelihood approach. 

For more details about the model, estimation, the finite sample behaviour in comparison to existing models (which ignore the presence of stayers or assume a fixed stayer fraction) and a read data application on Italian student mobility see Musta and Vittorietti (2025).

Here we demonstrate the method for generated data (under the setting 1 of the paper). The file 'main.R' contains the main estimation procedure, including construction of bootstrap confidence intervals. The file 'functions.R' contains the required functions for estimation, while the file "data-generation.R" is used to generate longitudinal data in the required format for the main function. 

References

E. Musta and M. Vittorietti (2025) A mover and stayer model with time-dependent stayer fraction. arXiv:2505.10065


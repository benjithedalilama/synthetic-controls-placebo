# synthetic-controls-placebo

Use this library to run a placebo test and replicate part of the paper “Estimating the population-level of vaccines using synthetic controls” to see if the claims in the paper hold. 

This uses the CausalImpact library at different intervention dates to see if there were consistent significant effects of PCVs, either positive or negative, in Brazil on all-cause pneumonia cases.

Just run the `synthetic.R` script and wait a while for the CausalImpact library to run. You will get graphs and impact objects which include absolute and relative effects along with confidence intervals in p-values.

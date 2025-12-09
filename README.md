# MathModeling_GlucoseInsulinBetaCells
Ordinary differential equations models written in Matlab describing dynamics of pancreatic beta-cells (immature &amp; mature populations) and plasma concentrations of glucose and insulin. Explore the effects of various regulatory control of beta-cells. Includes functions to generate plots for single-parameter sensitivity analysis.

# Description
Part of my PhD thesis work. My dissertation can be found here: https://escholarship.org/uc/item/8q2436np


# Usage Instructions
Run main.m. The variable model_structure is a string specifying the regulatory control of beta-cells. To explore other regulatory controls, model_structure string formatting details are found within the file get_regulatory_functions.m.
To use this code for your own system of ordinary differential equations (ODEs),
1. Create a new function file for your system of ODEs, based on ODE_GlucoseInsulin.m and ODE_RegulateBetaCells.m.
2. In the file get_parameters.m, add a new block of code similar to lines 60-90, "if strcmp(varargin{1}, "YourModelName")... end" that returns your model's parameters, etc. when calling get_parameters("YourModelName").
3. Depending on your model, especially if it has either 4 OR more than 6 dynamic variables, you may need to modify Y_transform.m. 


# Visuals 
Figure: Perturbing the value of the insulin resistance parameter. 
<img width="1380" height="1034" alt="Perturbing_KI" src="https://github.com/user-attachments/assets/32fc3d58-7f76-4c03-8aac-3573ad8e14b5" />


# Support/Help
Contact maggie.myers95@gmail.com


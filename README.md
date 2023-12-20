# stochastic-model-YAP-activation
A stochastic model for YAP activation within focal adhesion complexes

This repository contains the required files to run the 3D particle-based stochastic simulation for YAP activation within the focal adhesion complexes.

## Stochastic Simulation Algorithm (stochastic_model.m)
The stochastic_model.m file contains the simulation algorithm. The box features periodic boundary conditions in x-y directions, and the membrane is represented as a 2D plane at the bottom of the box. Initially, adhesions are randomly positioned on the membrane.

YAP molecules diffuse freely in the box. When they come into proximity with adhesions, they can bind to them, undergo phosphorylation, and subsequently detach from the adhesions. The simulation incorporates specific rates for reaction and diffusion events.

To visualize the simulation box, uncomment the line within the code where the plotPosition function is called.

## Plotting Results (plot_results.m)
The plot_results.m generates plots for pYAP ratio (phosphorylated YAP divided by total YAP) as a function of time. 


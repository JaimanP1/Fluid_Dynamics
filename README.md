# Nonlinear water waves research

Summer '23 research with Professor W. Choi, NJIT

Investigating dynamics of a contained fluid subject to a vertical oscillation while also including nonlinear contributions of surface tension and viscosity 

## How this works

The WaveAnimation.m file does exactly what the title suggests, produces a wave animation

The System folder contains files which work by passing a coupled nonlinear, second order system of differential equations through a Runge-Kutta routine

1. The system4implement_v4.m file gives control over the input, loop, and output for the routine
2. The system4_v6.m file contains the actual system
3. The rk4SingleStep.m file contains the algorithm for the fourth-order Runge-Kutta method
4. The l_vs_h_v2.m file plots physical quantities of interest
5. The ParametersClass.m file is a class which contains functions for all relevant variables, which can be called

### Envelope

The EnvelopeSystem folder contains files to plot an envelope for the solutions of the above system for given initial conditions. These files must be run first in order to generate solutions to the files in the System folder, namely system4implement_v4.m, which directly calls the results

1. The EnvelopeImplement_v4.m file gives control over the input, loop, and output for the routine
2. The EnvelopeSystem_v3.m contains the actual system
3. As above, the files work by calling the rk4 method and parameters from respective files
# HeatSolver.jl

A simple type-flexible model for heat diffusion through a soil column.

The model is taken from [1] and is of interest as it gives a good example of *stagnation* when integrated in single (rather than double) precision with round-to-nearest (heat does not diffuse effectively in single precision with the small choice of time-step). 

The effects of rounding error here can be mitigated here by either of (i) increasing the time-step Δt, or (ii) using a stochastic rounding scheme. 

Written in Julia, the code integrates the 1D heat equation dT/dt = D d^2T/dz^2 using a finite difference scheme (forward-Euler in time, centred finite difference in space).

The default parameters are as follows. 
Boundary conditions T(0,t) = 280 K at soil top, thermal insulation dT/dz(H,t) = 0 at soil bottom. 
Initial condition T(z,0) = 273.3 K. 
Thermal diffusivity D = 7 x 10^-7 m^2s^-1.
Soil depth H = 60 m.
Forecast length T = 60 years. 
Resolution Δt = 1800 s, Δz = 1 m.

# Usage

A simple usage example coming soon...

# References

[1] Dawson et. al - *Reliable low precision simulations in land surface models* (Climate Dynamics 51:2657–2666, 2018)
